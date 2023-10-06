import pytest

from quartic_solver import Solvers


@pytest.mark.parametrize(
    "input_coeffs, expected_roots",
    [
        ([-6, 0, 5, 0, 2], [-0.941, 0.941]),
    ],
)
def test_solve_biquadratic(input_coeffs, expected_roots):
    roots = Solvers.solve_quartic(input_coeffs)

    print(roots)
    for _, root in roots:
        assert any(
            pytest.approx(root, rel=1e-6) == expected_root
            for expected_root in expected_roots
        )


if __name__ == "__main__":
    test_solve_biquadratic([6, 0, 5, 0, -2], [-0.941, 0.941])
