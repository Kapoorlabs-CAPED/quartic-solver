import pytest

from quartic_solver import Solvers


# Coefficients are in the order of constant, x^1, x^2, x^3, x^4
@pytest.mark.parametrize(
    "input_coeffs, expected_roots",
    [
        ([6, -5, 1, 0, 0], [2, 3]),
        ([2, 0, -4, 3, 0], [-0.588911]),
        ([2, 0, -4, 0, 1], [-0.76537, 0.76537, 1.84776, -1.84776]),
        ([6, 5, -6, -3, 2], [1.5, 2, -1]),
    ],
)
def test_solve_quartic(input_coeffs, expected_roots):
    roots = Solvers.solve_quartic(input_coeffs)

    for _, root in roots:
        print(root)
        assert any(
            pytest.approx(root, rel=1e-4) == expected_root
            for expected_root in expected_roots
        )


if __name__ == "__main__":
    pytest.main()
