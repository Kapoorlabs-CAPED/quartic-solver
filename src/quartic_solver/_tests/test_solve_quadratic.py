import pytest

from quartic_solver import Solvers


@pytest.mark.parametrize(
    "input_coeffs, expected_roots",
    [
        ([6, -5, 1], [2, 3]),
    ],
)
def test_solve_quadratic(input_coeffs, expected_roots):
    roots = Solvers.solve_quadratic(input_coeffs)

    for _, root in roots:
        assert any(
            pytest.approx(root, rel=1e-6) == expected_root
            for expected_root in expected_roots
        )


if __name__ == "__main__":
    test_solve_quadratic([6, -5, 1], [2, 3])
