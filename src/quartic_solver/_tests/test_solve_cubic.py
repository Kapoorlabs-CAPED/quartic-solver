import pytest

from quartic_solver import Solvers


@pytest.mark.parametrize(
    "input_coeffs, expected_roots",
    [
        ([-6, 5, 0, 2], [0.935]),
    ],
)
def test_solve_cubic(input_coeffs, expected_roots):
    roots = Solvers.solve_cubic(input_coeffs)

    print(roots)
    for _, root in roots:
        assert any(
            pytest.approx(root, rel=1e-6) == expected_root
            for expected_root in expected_roots
        )


if __name__ == "__main__":
    test_solve_cubic([-6, 0, 5, 2], [0.935])
