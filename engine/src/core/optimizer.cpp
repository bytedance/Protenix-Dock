/*
* Copyright (C) 2025 ByteDance and/or its affiliates
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "bytedock/core/optimizer.h"

#include <cmath>
#include <deque>

#include "bytedock/core/interaction.h"
#include "bytedock/ext/logging.h"

namespace bytedock {

bool serial_optimizer::minimize_energy(binding_input& in, optimized_result& out) {
    out.ligand_xyz = in.ligand->get_pose(in.pose_id);  // Copy
    out.torsions.resize(in.receptor->num_torsions(), 0_r);
    binding_system_interactions model(in.receptor, in.ligand, in.cache);
    bool converged = step_.apply(model, out);
    out.receptor_xyz = in.receptor->apply_parameters(out.torsions);  // Move
    out.offset = in.pose_id;
    return converged;
}

void serial_optimizer::fill(blocking_queue<name_and_task>& out_queue) {
    name_and_task pair = std::move(task_queue_.pop());
    while (!task_queue_.is_eoq(pair)) {
        LOG_INFO << "Task [" << pair.first << "] has started...";
        {
            step_timer t(REFINE_BINDING_POSE_STEP);
            if (!minimize_energy(pair.second.feed, pair.second.fetch)) {
                LOG_WARNING << "Minimization of [" << pair.first
                            << "] failed to converge!";
            }
        }
        out_queue.push(std::move(pair));
        pair = std::move(task_queue_.pop());
    }
    task_queue_.close();  // Put back the EOF flag
    out_queue.close();
}

inline std::vector<param_t> flat_gradients(
    const std::vector<param_t>& torsion_gradient,
    const molecule_pose& ligand_gradient
) {
    std::vector<param_t> merged;
    merged.resize(torsion_gradient.size() + ligand_gradient.size() * 3);
    size_t offset = 0;
    for (const auto& item: ligand_gradient) {
        merged[offset] = item.xyz[0];
        offset++;
        merged[offset] = item.xyz[1];
        offset++;
        merged[offset] = item.xyz[2];
        offset++;
    }
    for (const auto& item : torsion_gradient) {
        merged[offset] = item;
        offset++;
    }
    return merged;
}

inline param_t get_abs_max(const std::vector<param_t>& array) {
    param_t value = 0., tmp;
    for (auto item : array) {
        tmp = std::abs(item);
        if (tmp > value) value = tmp;
    }
    return value;
}

inline param_t get_abs_sum(const std::vector<param_t>& array) {
    param_t value = 0.;
    for (auto item : array) value += std::abs(item);
    return value;
}

inline std::vector<param_t> get_negative(const std::vector<param_t>& vin) {
    std::vector<param_t> vout(vin.size());
    for (size_t i = 0; i < vout.size(); i++) vout[i] = -vin[i];
    return vout;
}

inline std::vector<param_t> minus(const std::vector<param_t>& left,
                                  const std::vector<param_t>& right) {
    std::vector<param_t> vout(left.size());
    for (size_t i = 0; i < right.size(); i++) vout[i] = left[i] - right[i];
    return vout;
}

inline std::vector<param_t> scale(const std::vector<param_t>& vin,
                                  const param_t factor) {
    std::vector<param_t> vout(vin.size());
    for (size_t i = 0; i < vout.size(); i++) vout[i] = vin[i] * factor;
    return vout;
}

inline param_t dot_product(const std::vector<param_t>& left,
                           const std::vector<param_t>& right) {
    param_t value = 0.;
    for (size_t i = 0; i < right.size(); i++) value += left[i] * right[i];
    return value;
}

inline void add_inplace(std::vector<param_t>& array,
                        const std::vector<param_t>& other,
                        const param_t factor) {
    for (size_t i = 0; i < array.size(); i++) array[i] += other[i] * factor;
}

inline void add_grad(const param_t t, const std::vector<param_t>& d,
                     std::vector<param_t>& torsions, molecule_pose& ligand_xyz) {
    size_t offset = 0;
    for (size_t i = 0; i < ligand_xyz.size(); i++) {
        ligand_xyz[i].xyz[0] += d[offset] * t;
        offset++;
        ligand_xyz[i].xyz[1] += d[offset] * t;
        offset++;
        ligand_xyz[i].xyz[2] += d[offset] * t;
        offset++;
    }
    for (size_t i = 0; i < torsions.size(); i++) {
        torsions[i] += d[offset] * t;
        offset++;
    }
}

class obj_func_closure {
public:
    explicit obj_func_closure(
        binding_system_interactions& model, size_t& fun_evals,
        std::vector<param_t>& receptor_torsions, molecule_pose& ligand_xyz,
        std::vector<param_t>& torsion_gradient, molecule_pose& ligand_gradient
    ) : bsi_(model), evals_(fun_evals), rt_(receptor_torsions), lxyz_(ligand_xyz),
        tg_(torsion_gradient), lg_(ligand_gradient) {
        init_rt_ = rt_;  // Copy
        init_lxyz_ = lxyz_;  // Copy
    }

    param_t evaluate_direction(const param_t t,
                               const std::vector<param_t>& d,
                               std::vector<param_t>& flat_grad) {
        add_grad(t, d, rt_, lxyz_);
        bsi_.reset_gradients(tg_, lg_);
        param_t energy = bsi_.put_gradients(rt_, lxyz_, tg_, lg_);
        flat_grad = flat_gradients(tg_, lg_);  // Move
        ++evals_;
        rt_ = init_rt_;
        lxyz_ = init_lxyz_;
        bsi_.reset_gradients(tg_, lg_);
        return energy;
    }

    std::string format_parameters(const std::string prefix = "") {
        std::ostringstream oss;
        oss << "{" << std::endl;
        oss << prefix << "  \"ligand_xyz\": [" << std::endl;
        const param_t* ptr;
        for (size_t i = 0; i < lxyz_.size() - 1; i++) {
            ptr = lxyz_[i].xyz;
            oss << prefix << "    [" << ptr[0] << ", " << ptr[1] << ", " << ptr[2] << "]," << std::endl;
        }
        ptr = lxyz_[lxyz_.size() - 1].xyz;
        oss << prefix << "    [" << ptr[0] << ", " << ptr[1] << ", " << ptr[2] << "]" << std::endl;
        oss << prefix << "]," << std::endl;
        oss << prefix << "  \"torsions\": [";
        for (size_t  i = 0; i < rt_.size() - 1; i++) {
            oss << rt_[i] << ", ";
        }
        oss << rt_[rt_.size() - 1];
        oss << prefix << "  ]" << std::endl;
        oss << prefix << "}";
        return oss.str();
    }

private:
    binding_system_interactions& bsi_;
    size_t& evals_;

    std::vector<param_t>& rt_;
    molecule_pose& lxyz_;

    std::vector<param_t> init_rt_;
    molecule_pose init_lxyz_;

    std::vector<param_t>& tg_;
    molecule_pose& lg_;
};

inline param_t cubic_interpolate(const param_t x1, const param_t f1, const param_t g1,
                                 const param_t x2, const param_t f2, const param_t g2,
                                 const param_t xmin_bound, const param_t xmax_bound) {
    /**
     * Code for most common case: cubic interpolation of 2 points w/ function and
     * derivative values for both Solution in this case (where x2 is the farthest
     * point):
     * d1 = g1 + g2 - 3_r * (f1 - f2) / (x1 - x2)
     * d2 = sqrt(d1^2 - g1 * g2)
     * min_pos = x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2_r * d2))
     * t_new = min(max(min_pos, xmin_bound), xmax_bound)
     */
    param_t d1 = g1 + g2 - 3_r * (f1 - f2) / (x1 - x2);
    param_t d2_square = d1 * d1 - g1 * g2;
    param_t d2;
    if (d2_square >= 0_r) {
        d2 = std::sqrt(d2_square);
        // `d1` => `min_pos`
        if (x1 <= x2) {
            d1 = x2 - (x2 - x1) * ((g2 + d2 - d1) / (g2 - g1 + 2_r * d2));
        } else {
            d1 = x1 - (x1 - x2) * ((g1 + d2 - d1) / (g1 - g2 + 2_r * d2));
        }
        return std::min(MATH_PAIR_MAX(d1, xmin_bound), xmax_bound);
    } else {
        return (xmin_bound + xmax_bound) / 2_r;
    }

}

static void strong_wolfe(obj_func_closure& obj_func,
                         param_t& t,  // Both in & out
                         const std::vector<param_t>& d,
                         param_t& f,  // Both in & out
                         std::vector<param_t>& g,  // Both in & out
                         const param_t gtd,
                         const param_t gtd_threshold = 1e5,
                         const param_t c1 = 1e-4,
                         const param_t c2 = 0.9,
                         const param_t tolerance_change = 1e-9,
                         const size_t max_ls = 25) {
    param_t bracket[2], bracket_f[2], bracket_gtd[2];
    std::vector<param_t> bracket_g[2];

    // Evaluate objective and gradient using initial step
    std::vector<param_t> g_new;
    param_t f_new = obj_func.evaluate_direction(t, d, g_new);
    if (std::isnan(f_new)) {
        f = f_new;
        return;
    }
    param_t gtd_new = dot_product(g_new, d);

    // Wolfe condition along large descent direction not satisfied
    if (std::abs(gtd_new) > gtd_threshold) {
        f = f_new;
        g = std::move(g_new);
        return;
    }

    // Bracket an interval containing a point satisfying the Wolfe criteria
    param_t t_prev = 0., f_prev = f, gtd_prev = gtd, min_step, max_step, tmp;
    std::vector<param_t> g_prev = g;  // Copy
    bool done = false;
    size_t ls_iter = 0;
    while (ls_iter < max_ls) {
        LOG_DEBUG << "strong_wolfe/normal_case/ls_iter = " << ls_iter;
        if ((f_new > f + c1 * t * gtd) || (ls_iter > 1 && f_new >= f_prev)) {
            bracket[0] = t_prev;
            bracket[1] = t;
            bracket_f[0] = f_prev;
            bracket_f[1] = f_new;
            bracket_g[0] = std::move(g_prev);
            bracket_g[1] = std::move(g_new);
            bracket_gtd[0] = gtd_prev;
            bracket_gtd[1] = gtd_new;
            break;
        }
        if (std::abs(gtd_new) <= -c2 * gtd) {
            bracket[0] = t;
            bracket_f[0] = f_new;
            bracket_g[0] = std::move(g_new);
            done = true;
            break;
        }
        if (gtd_new >= 0) {
            bracket[0] = t_prev;
            bracket[1] = t;
            bracket_f[0] = f_prev;
            bracket_f[1] = f_new;
            bracket_g[0] = std::move(g_prev);
            bracket_g[1] = std::move(g_new);
            bracket_gtd[0] = gtd_prev;
            bracket_gtd[1] = gtd_new;
            break;
        }

        // Interpolate
        min_step = t + 0.01_r * (t - t_prev);
        max_step = t * 10_r;
        tmp = t;
        t = cubic_interpolate(t_prev, f_prev, gtd_prev, t, f_new, gtd_new,
                              min_step, max_step);

        // Next step
        t_prev = tmp;
        f_prev = f_new;
        g_prev = std::move(g_new);
        gtd_prev = gtd_new;
        f_new = obj_func.evaluate_direction(t, d, g_new);
        gtd_new = dot_product(g_new, d);
        ls_iter++;
    }

    // Reached max number of iterations?
    param_t d_norm = 0.;
    if (ls_iter == max_ls) {
        bracket[0] = 0;
        bracket[1] = t;
        bracket_f[0] = f;
        bracket_f[1] = f_new;
        bracket_g[0] = std::move(g);
        bracket_g[1] = std::move(g_new);
    } else {
        d_norm = get_abs_max(d);
    }

    /**
     * Zoom phase: we now have a point satisfying the criteria, or a bracket around it.
     * We refine the bracket until we find the exact point satisfying the criteria.
     */
    param_t eps;
    bool insuf_progress = false;
    size_t low_pos = 0, high_pos = 1;
    if (!done && bracket_f[0] > bracket_f[1]) {
        low_pos = 1;
        high_pos = 0;
    }
    while (!done && ls_iter < max_ls) {
        LOG_DEBUG << "strong_wolfe/zoom_phase/ls_iter = " << ls_iter;
        if (bracket[0] > bracket[1]) {
            min_step = bracket[1];
            max_step = bracket[0];
        } else {
            min_step = bracket[0];
            max_step = bracket[1];
        }

        // Line-search bracket is so small
        eps = max_step - min_step;
        if (eps * d_norm < tolerance_change) break;

        // Compute new trial value
        t = cubic_interpolate(bracket[0], bracket_f[0], bracket_gtd[0],
                              bracket[1], bracket_f[1], bracket_gtd[1],
                              min_step, max_step);

        /**
         * Test that we are making sufficient progress. In case `t` is so close to
         * boundary, we mark that we are making insufficient progress, and if we have
         * made insufficient progress in the last step, or `t` is at one of the
         * boundary, we will move `t` to a position which is `0.1 * len(bracket)` away
         * from the nearest boundary point.
         */
        eps *= 0.1_r;
        if (std::min(max_step - t, t - min_step) < eps) {
            if (insuf_progress || t >= max_step || t <= min_step) {
                if (std::abs(t - max_step) < std::abs(t - min_step)) {
                    t = max_step - eps;
                } else {
                    t = min_step + eps;
                }
                insuf_progress = false;
            } else {
                insuf_progress = true;
            }
        } else {
            insuf_progress = false;
        }

        // Evaluate new point
        f_new = obj_func.evaluate_direction(t, d, g_new);
        gtd_new = dot_product(g_new, d);
        ls_iter++;

        if (f_new > (f + c1 * t * gtd) || f_new >= bracket_f[low_pos]) {
            // Armijo condition not satisfied or not lower than lowest point
            bracket[high_pos] = t;
            bracket_f[high_pos] = f_new;
            bracket_g[high_pos] = std::move(g_new);
            bracket_gtd[high_pos] = gtd_new;
            if (bracket_f[0] > bracket_f[1]) {
                low_pos = 1;
                high_pos = 0;
            } else {
                low_pos = 0;
                high_pos = 1;
            }
        } else {
            if (std::abs(gtd_new) <= -c2 * gtd) {
                // Wolfe conditions satisfied
                done = true;
            } else if (gtd_new * (bracket[high_pos] - bracket[low_pos]) >= 0_r) {
                // Old high becomes new low
                bracket[high_pos] = bracket[low_pos];
                bracket_f[high_pos] = bracket_f[low_pos];
                bracket_g[high_pos] = std::move(bracket_g[low_pos]);
                bracket_gtd[high_pos] = bracket_gtd[low_pos];
            }

            // New point becomes new low
            bracket[low_pos] = t;
            bracket_f[low_pos] = f_new;
            bracket_g[low_pos] = std::move(g_new);
            bracket_gtd[low_pos] = gtd_new;
        }
    }

    // Return stuff
    t = bracket[low_pos];
    f = bracket_f[low_pos];
    g = std::move(bracket_g[low_pos]);
}

bool lbfgs_step::apply(binding_system_interactions& model, optimized_result& out) {
    std::vector<param_t> torsion_gradient;
    molecule_pose ligand_gradient;
    model.init_gradients(torsion_gradient, ligand_gradient);

    // Evaluate initial f(x) and df/dx
    out.energy = model.put_gradients(out.torsions, out.ligand_xyz,
                                     torsion_gradient, ligand_gradient);
    out.nevals = 1;

    // Meet optimal condition
    std::vector<param_t> flat_grad = flat_gradients(torsion_gradient, ligand_gradient);
    if (get_abs_max(flat_grad) < tolerance_grad_) return true;

    // Decalre variables for tracing
    bool high_grad_change = false, converged = false, updated = false;
    std::vector<param_t> prev_flat_grad, y, d, s, q;
    param_t loss = out.energy, prev_loss;
    param_t H_diag, t, ys, be_i, gtd;
    std::deque<std::vector<param_t> > old_dirs, old_stps;
    std::deque<param_t> ro;
    std::vector<param_t> al(history_size_);
    size_t n_iter = 0, num_old;

    // Optimize for a max of `max_iter_` iterations
    while (n_iter < max_iter_ && out.nevals < max_eval_) {
        LOG_DEBUG << "l-bfgs-gl/n_iter = " << n_iter;
        n_iter++;  // keep track of nb of iterations

        // --- Compute gradient descent direction -->
        if (n_iter > 1) {
            y = minus(flat_grad, prev_flat_grad);
            high_grad_change = get_abs_max(y) > 1e8_r;
        }
        if (n_iter == 1 || high_grad_change) {
            d = get_negative(flat_grad);  // Move assignment
            old_dirs.clear();
            old_stps.clear();
            ro.clear();
            H_diag = 1;
        } else {
            // Do L-BFGS update (along with memory)
            s = scale(d, t);
            ys = dot_product(y, s);
            if (ys > 1e-10_r) {  // Update memory
                if (old_dirs.size() == history_size_) {
                    old_dirs.pop_front();
                    old_stps.pop_front();
                    ro.pop_front();
                }

                // Update scale of initial Hessian approximation before `y` is moved
                H_diag = ys / dot_product(y, y);

                // Store new direction/step
                old_dirs.push_back(std::move(y));
                old_stps.push_back(std::move(s));
                ro.push_back(1_r / ys);
            }

            /**
             * Compute the approximate (L-BFGS) inverse Hessian multiplied by the
             * gradient.
             */
            num_old = old_dirs.size();

            // # iteration in L-BFGS loop collapsed to use just one buffer
            q = get_negative(flat_grad);
            for (int i = num_old - 1; i > -1; i--) {
                al[i] = dot_product(old_stps[i], q) * ro[i];
                add_inplace(q, old_dirs[i], -al[i]);
            }

            // Multiply by initial Hessian to get the final direction
            d = scale(q, H_diag);
            for (size_t i = 0; i < num_old; i++) {
                be_i = dot_product(old_dirs[i], d) * ro[i];
                add_inplace(d, old_stps[i], al[i] - be_i);
            }
        }

        prev_flat_grad = flat_grad;  // Copy
        prev_loss = loss;

        // --- Compute step length -->
        if (n_iter == 1 or high_grad_change) {  // Reset initial guess for step size
            t = std::min(1_r, 1_r / get_abs_sum(flat_grad)) * lr_ * 0.1_r;
        } else {
            t = lr_;
        }

        // Directional derivative
        gtd = dot_product(flat_grad, d);
        if (gtd > -tolerance_change_) {
            converged = true;
            break;
        }

        // Perform line search
        obj_func_closure closure(model, out.nevals,
                                 out.torsions, out.ligand_xyz,
                                 torsion_gradient, ligand_gradient);
        strong_wolfe(closure, t, d, loss, flat_grad, gtd);
        if (std::isnan(loss)) break;  // Not converged but invalid
        add_grad(t, d, out.torsions, out.ligand_xyz);
        updated = true;
        LOG_DEBUG << "l-bfgs-gl/current_evals = " << out.nevals;

        // --- Check conditions -->
        if (get_abs_max(flat_grad) <= tolerance_grad_) {
            converged = true;
            break;
        }
        s = scale(d, t);
        if (get_abs_max(s) <= tolerance_change_) {
            converged = true;
            break;
        }
        if (std::abs(loss - prev_loss) < tolerance_change_) {
            converged = true;
            break;
        }
    }

    // Generate outputs
    if (updated) {
        out.energy = model.put_gradients(out.torsions, out.ligand_xyz,
                                         torsion_gradient, ligand_gradient);
        ++out.nevals;
    }
    return converged;
}

}
