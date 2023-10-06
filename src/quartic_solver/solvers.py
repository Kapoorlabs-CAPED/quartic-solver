import math
from typing import List, Tuple


class Solvers:
    @staticmethod
    def solve_quartic(p: List[float]) -> List[Tuple[int, float]]:
        root_map = []

        rat2 = 2
        rat3 = 3
        rat4 = 4
        rat6 = 6

        q0 = p[0] / p[4]
        q1 = p[1] / p[4]
        q2 = p[2] / p[4]
        q3 = p[3] / p[4]

        q3_fourth = q3 / rat4
        q3_fourth_sqr = q3_fourth * q3_fourth
        c0 = q0 - q3_fourth * (q1 - q3_fourth * (q2 - q3_fourth_sqr * rat3))
        c1 = q1 - rat2 * q3_fourth * (q2 - rat4 * q3_fourth_sqr)
        c2 = q2 - rat6 * q3_fourth_sqr

        root_map_local = Solvers.solve_depressed_quartic(c0, c1, c2)
        for rm in root_map_local:
            root = rm[1] - q3_fourth
            root_map.append((rm[0], root))

        return root_map

    @staticmethod
    def solve_depressed_quartic(
        c0: float, c1: float, c2: float
    ) -> List[Tuple[int, float]]:
        root_map = []

        zero = 0.0
        if c0 == zero:
            root_map_local = Solvers.solve_depressed_cubic(c1, c2)

            for rm in root_map_local:
                if rm[1] != zero:
                    # The cubic does not have a root of zero, Insert one for quartic
                    root_map.append((1, zero))

            return root_map

        if c1 == zero:
            root_map = Solvers.solve_biquadratic(c0, c2)
            return root_map

        rat2, rat3, rat4, rat8, rat12, rat16, rat27, rat36 = (
            2.0,
            3.0,
            4.0,
            8.0,
            12.0,
            16.0,
            27.0,
            36.0,
        )
        c0sqr, c1sqr, c2sqr = c0 * c0, c1 * c1, c2 * c2
        delta = c1sqr * (
            -rat27 * c1sqr + rat4 * c2 * (rat36 * c0 - c2sqr)
        ) + rat16 * c0 * (c2sqr * (c2sqr - rat8 * c0) + rat16 * c0sqr)
        a0, a1 = rat12 * c0 + c2sqr, rat4 * c0 - c2sqr

        if delta > zero:
            print("Delta > zero QuarticDep", delta)
            if c2 < zero and a1 < zero:
                # Four simple roots
                root_map_local = Solvers.solve_cubic(
                    [c1sqr - rat4 * c0 * c2, rat8 * c0, rat4 * c2, -rat8]
                )
                t = root_map_local[-1][1]
                alpha_sqr = rat2 * t - c2
                alpha = alpha_sqr**0.5
                sgn_c1 = 1.0 if c1 > zero else -1.0
                arg = t * t - c0
                beta = sgn_c1 * (arg**0.5)
                D0 = alpha_sqr - rat4 * (t + beta)
                sqrtD0 = max(D0, zero) ** 0.5
                D1 = alpha_sqr - rat4 * (t - beta)
                sqrtD1 = max(D1, zero) ** 0.5
                root0 = (alpha - sqrtD0) / rat2
                root1 = (alpha + sqrtD0) / rat2
                root2 = (-alpha - sqrtD1) / rat2
                root3 = (-alpha + sqrtD1) / rat2
                root_map.append((1, root0))
                root_map.append((1, root1))
                root_map.append((1, root2))
                root_map.append((1, root3))
            else:
                # c2 >= 0 or a1 >= 0, Roots are complex
                pass
            return root_map
        elif delta < zero:
            print("Delta < zero QuarticDep")
            # Two simple real roots, one complex conjugate pair
            root_map_local = Solvers.solve_cubic(
                [c1sqr - rat4 * c0 * c2, rat8 * c0, rat4 * c2, -rat8]
            )
            t = root_map_local[-1][1]
            alpha_sqr = rat2 * t - c2
            alpha = alpha_sqr**0.5
            sgn_c1 = 1.0 if c1 > zero else -1.0
            arg = t * t - c0
            beta = sgn_c1 * (arg**0.5)
            if sgn_c1 > 0:
                print("sgn_c1 > 0 QuarticDep")
                D1 = alpha_sqr - rat4 * (t - beta)
                sqrtD1 = max(D1, zero) ** 0.5
                root0 = (-alpha - sqrtD1) / rat2
                root1 = (-alpha + sqrtD1) / rat2
            else:
                print("sgn_c1 < 0 QuarticDep")
                D0 = alpha_sqr - rat4 * (t + beta)
                sqrtD0 = max(D0, zero) ** 0.5
                root0 = (alpha - sqrtD0) / rat2
                root1 = (alpha + sqrtD0) / rat2
            root_map.append((1, root0))
            root_map.append((1, root1))
            return root_map
        else:  # delta == 0
            if a1 > zero or (c2 > zero and (a1 != zero or c1 != zero)):
                # One double real root
                rat9 = 9
                root0 = -c1 * a0 / (rat9 * c1sqr - rat2 * c2 * a1)
                root_map.append((2, root0))
            else:
                if a0 != zero:
                    # One double real root, two simple real roots
                    rat9 = 9
                    root0 = -c1 * a0 / (rat9 * c1sqr - rat2 * c2 * a1)
                    alpha = rat2 * root0
                    beta = c2 + rat3 * root0 * root0
                    discr = alpha * alpha - rat4 * beta
                    temp1 = max(discr, zero) ** 0.5
                    root1 = (-alpha - temp1) / rat2
                    root2 = (-alpha + temp1) / rat2
                    root_map.append((2, root0))
                    root_map.append((1, root1))
                    root_map.append((1, root2))
                else:
                    # One triple real root, one simple real root
                    root0 = -rat3 * c1 / (rat4 * c2)
                    root1 = -rat3 * root0
                    root_map.append((3, root0))
                    root_map.append((1, root1))
            return root_map

    @staticmethod
    def solve_depressed_quadratic(c0: float) -> List[Tuple[int, float]]:
        root_map = []

        zero = 0.0
        if c0 < zero:
            # Two simple roots
            root1 = (-c0) ** 0.5
            root0 = -root1
            root_map.append((1, root0))
            root_map.append((1, root1))
        elif c0 == zero:
            # One double root
            root_map.append((2, zero))
        else:
            # Roots are complex (not adding any roots in this case)
            pass

        return root_map

    @staticmethod
    def solve_depressed_cubic(c0: float, c1: float) -> List[Tuple[int, float]]:
        root_map = []

        zero = 0.0
        if c0 == zero:
            root_map_local = Solvers.solve_depressed_quadratic(c1)

            for rm in root_map_local:
                if rm[1] != 0.0:
                    # The quadratic does not have a root of zero, Insert one for cubic
                    root_map.append((1, zero))

            return root_map

        one_third = 1.0 / 3.0

        if c1 == zero:
            if c0 > zero:
                root0 = -(c0**one_third)
            else:
                root0 = (-c0) ** one_third

            root_map.append((1, root0))
            return root_map

        rat2 = 2.0
        rat3 = 3.0
        rat4 = 4.0
        rat27 = 27.0
        rat108 = 108.0

        delta = -(rat4 * c1 * c1 * c1 + rat27 * c0 * c0)
        if delta > zero:
            print("Delta > zero")
            # Three simple roots
            delta_div_108 = delta / rat108
            beta_re = -c0 / rat2
            beta_im = (delta_div_108) ** 0.5
            theta = math.atan2(beta_im, beta_re)
            theta_div_3 = theta / rat3
            angle = theta_div_3
            cs = math.cos(angle)
            sn = math.sin(angle)
            rho_sqr = beta_re * beta_re + beta_im * beta_im
            rho_pow_third = rho_sqr ** (1.0 / 6.0)
            temp0 = rho_pow_third * cs
            temp1 = rho_pow_third * sn * (3**0.5)
            root0 = rat2 * temp0
            root1 = -temp0 - temp1
            root2 = -temp0 + temp1
            root_map.append((1, root0))
            root_map.append((1, root1))
            root_map.append((1, root2))
        elif delta < zero:
            print("Delta < zero")
            # One Simple root
            delta_div_108 = delta / rat108
            temp0 = -c0 / rat2
            temp1 = (-delta_div_108) ** 0.5
            temp2 = temp0 - temp1
            temp3 = temp0 + temp1

            if temp2 >= zero:
                temp22 = temp2**one_third
            else:
                temp22 = (-temp2) ** one_third

            if temp3 >= zero:
                temp33 = temp3**one_third
            else:
                temp33 = (-temp3) ** one_third

            root0 = temp22 + temp33
            root_map.append((1, root0))
        else:
            print("Delta Nothing zero")
            # One simple root and one double root.
            root0 = -rat3 * c0 / (rat2 * c1)
            root1 = -rat2 * root0
            root_map.append((2, root0))
            root_map.append((1, root1))

        return root_map

    @staticmethod
    def solve_biquadratic(c0, c2):
        Rootmap = []

        zero = 0
        rat2 = 2
        rat256 = 256
        c2_half = c2 / rat2
        a1 = c0 - c2_half * c2_half
        delta = rat256 * c0 * a1 * a1
        if delta > zero:
            if c2 < zero:
                if a1 < zero:
                    # Four simple roots
                    temp0 = (-a1) ** 0.5
                    temp1 = -c2_half - temp0
                    temp2 = -c2_half + temp0
                    root0 = (temp1) ** 0.5
                    root1 = -root0
                    root2 = (temp2) ** 0.5
                    root3 = -root2
                    Rootmap.append((1, root0))
                    Rootmap.append((1, root1))
                    Rootmap.append((1, root2))
                    Rootmap.append((1, root3))
                else:
                    # Roots are complex
                    pass
            else:
                # c2 > 0
                # Roots are complex
                pass
        elif delta < zero:
            # Two simple real roots
            root0 = (-c2_half) ** 0.5
            root1 = -root0
            Rootmap.append((1, root0))
            Rootmap.append((1, root1))
        else:
            if c2 < zero:
                # Two double real roots
                root0 = (-c2_half) ** 0.5
                root1 = -root0
                Rootmap.append((2, root0))
                Rootmap.append((2, root1))
            else:
                # Roots are complex
                pass

        return Rootmap

    @staticmethod
    def solve_cubic(p):
        Rootmap = []

        rat2 = 2
        rat3 = 3

        q0 = p[0] / p[3]
        q1 = p[1] / p[3]
        q2 = p[2] / p[3]

        q2third = q2 / rat3
        c0 = q0 - q2third * (q1 - rat2 * q2third * q2third)
        c1 = q1 - rat2 * q2third * q2third
        RootmapLocal = Solvers.solve_depressed_cubic(c0, c1)

        for rm in RootmapLocal:
            root = rm[1] - q2third
            Rootmap.append((rm[0], root))

        return Rootmap

    @staticmethod
    def solve_quadratic(p):
        Rootmap = []

        q0 = p[0] / p[2]
        q1 = p[1] / p[2]

        disc = q1 * q1 - 4 * q0
        rootA = (-q1 + disc) / 2.0
        rootB = (-q1 - disc) / 2.0
        Rootmap.append((1, rootA))
        Rootmap.append((1, rootB))

        return Rootmap