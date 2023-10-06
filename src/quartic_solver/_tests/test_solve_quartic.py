import pytest

from quartic_solver import Solvers


@pytest.mark.parametrize(
    "input_coeffs, expected_roots",
    [
        ([2, -3, -6, 5, 6], [-0.716, 2.251]),
    ],
)
def test_solve_quartic(input_coeffs, expected_roots):
    roots = Solvers.solve_quartic(input_coeffs)

    for _, root in roots:
        print(root)
        assert any(
            pytest.approx(root, rel=1e-6) == expected_root
            for expected_root in expected_roots
        )


if __name__ == "__main__":
    test_solve_quartic([2, -3, -6, 5, 6], [-0.716, 2.251])
