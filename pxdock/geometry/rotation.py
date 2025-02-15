# Copyright (C) 2025 ByteDance and/or its affiliates

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import torch


def euler_rotation_matrix(
    angles: torch.Tensor,
):
    """
    Calculate the rotation matrix from Euler angles.
    For the theory of Euler angles, see: https://en.wikipedia.org/wiki/Euler_angles

    Args:
        angles: torch.tensor
            - Euler angles: alpha, beta, gamma
            - shape=[..., 3]

    Returns:
        R: torch.tensor
            - rotation matrix
            - shape=[..., 3, 3]
    """

    assert angles.shape[-1] == 3

    sinx = torch.sin(angles)
    cosx = torch.cos(angles)

    sa, sb, sg = sinx[..., 0], sinx[..., 1], sinx[..., 2]
    ca, cb, cg = cosx[..., 0], cosx[..., 1], cosx[..., 2]

    R = torch.stack(
        [
            ca * cb,
            ca * sb * sg - sa * cg,
            ca * sb * cg + sa * sg,
            sa * cb,
            sa * sb * sg + ca * cg,
            sa * sb * cg - ca * sg,
            -sb,
            cb * sg,
            cb * cg,
        ],
        dim=-1,
    )  # shape=[..., 9]

    # reshape last dim [9] --> [3, 3]
    R = R.view(*R.shape[:-1], 3, 3)

    return R


def rodrigues_rotation_matrix(u: torch.Tensor, theta: torch.Tensor):
    """
    Calculate the rotation matrix from Rodrigues formula.
    For the theory of Redrigues formula, refer to https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

    Args:
        u: torsion shaft, shape=[..., 3]
        theta: torch.tensor, shape [..., 1] or [...]

    Returns:
        R_matrix: shape [..., 3, 3]

    """

    assert u.shape[-1] == 3
    x, y, z = u[..., 0], u[..., 1], u[..., 2]

    if len(theta.shape) > len(x.shape):
        theta = theta.squeeze(dim=-1)
    assert theta.shape == x.shape
    st = torch.sin(theta)
    ct = torch.cos(theta)

    R = torch.stack(
        [
            ct + torch.pow(x, 2) * (1 - ct),
            x * y * (1 - ct) - z * st,
            y * st + x * z * (1 - ct),
            z * st + x * y * (1 - ct),
            ct + torch.pow(y, 2) * (1 - ct),
            -x * st + y * z * (1 - ct),
            -y * st + x * z * (1 - ct),
            x * st + y * z * (1 - ct),
            ct + torch.pow(z, 2) * (1 - ct),
        ],
        dim=-1,
    )
    R = R.view(*R.shape[:-1], 3, 3)

    return R
