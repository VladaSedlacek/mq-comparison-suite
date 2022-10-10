from pathlib import Path


def elt_to_str(q, a):
    if q.is_prime():
        # this is just a special case of the latter, but it is more readable
        return str(hex(a))[2:].zfill(2)
    else:
        d = log(q, radical(q))
        return str(hex(sum([2**i * a.polynomial()[i].lift() for i in range(d)])))[2:].zfill(2)


def str_to_elt(q, s, field=None):
    if field is not None:
        return field.fetch_int(int(s, 16))
    else:
        return GF(q).fetch_int(int(s, 16))


def delete_powers(eq):
    return sum([radical(mon) for mon in eq.monomials()])


def linear_combination(coefficients, objects):
    assert len(coefficients) == len(objects)
    return sum(c * o for c, o in zip(coefficients, objects))


def weil_decomposition(poly):
    if poly == 0:
        return []
    # Constant coefficients come first
    extension_coeffs = [c.polynomial().list() for c in poly.coefficients()]
    max_len = max(len(ec) for ec in extension_coeffs)
    # Pad the coefficient lists
    for ec in extension_coeffs:
        for _ in range(max_len - len(ec)):
            ec.append(0)
    base_coeff_list = zip(*extension_coeffs)
    base_polys = [linear_combination(coeffs, poly.monomials()) for coeffs in base_coeff_list]
    # Convert the base polys into actual polynomials over the base field
    ext_ring = poly.parent()
    base_field = ext_ring.base().base()
    base_ring = ext_ring.change_ring(base_field)
    return [base_ring(base_poly) for base_poly in base_polys]


def compute_degrevlex_mons(var_list):
    var_list = var_list + [1]
    degrevlex_mons = []
    # degrevlex corresponds to reading the upper part of the quadratic form matrix column-wise
    for j in range(len(var_list)):
        for i in range(j + 1):
            degrevlex_mons.append(var_list[i] * var_list[j])
    return degrevlex_mons


def load_fukuoka(fukuoka_path):
    # the input file format is the one used by https://www.mqchallenge.org/ (and coincides with cb_gpu)
    with open(fukuoka_path, 'r') as f:
        lines = f.readlines()

        # perform basic sanity checks on the input file
        assert "Galois Field : " in lines[0]
        assert "Number of variables (n) : " in lines[1]
        assert "Number of polynomials (m) : " in lines[2]
        assert "Seed : " in lines[3]
        assert "Order : " in lines[4]
        assert lines[5].strip() == ""
        assert lines[6].strip() == "*********************"

        # parse the header
        field_str = lines[0].split("Galois Field : ")[1].strip().split(" / ")
        if len(field_str) > 1:
            base_str = field_str[0].split("[x]")[0].strip()
            modulus_str = field_str[1].strip()
            ring = sage_eval(f"PolynomialRing({base_str}, 'x')")
            q = ring.characteristic()
            modulus = ring(modulus_str)
            d = modulus.degree()
            field = GF(q ^ d, name='z', modulus=modulus)
        else:
            field = sage_eval(field_str[0].strip())

        n = int(lines[1].split("Number of variables (n) : ")[1].strip())
        m = int(lines[2].split("Number of polynomials (m) : ")[1].strip())
        seed = int(lines[3].split("Seed : ")[1].strip())
        order = lines[4].split("Order : ")[1].strip()

        assert order == 'graded reverse lex order'
        # other orders are not supported yet

        # construct the equations
        R = PolynomialRing(field, ['x%s' % p for p in range(1, n + 1)], order='degrevlex')
        degrevlex_mons = compute_degrevlex_mons(list(R.gens()))
        equations = []
        for line in lines[7:]:
            assert line.strip()[-1] == ';'
            str_coeffs = line.strip().split(" ")[:-1]
            coeffs = [str_to_elt(field.order, str_coeff, field) for str_coeff in str_coeffs]
            equations.append(linear_combination(coeffs, degrevlex_mons))
        assert len(equations) == m
        return EquationSystem(equations, seed=seed)


def convert_fukuoka_to_others(fukuoka_path):
    folder = Path(fukuoka_path).parent
    base_system_name = Path(fukuoka_path).stem
    FukuokaSystem = load_fukuoka(fukuoka_path)
    FukuokaSystem.save_all(folder, base_system_name)
    # for systems over GF(2), the new ".cb_gpu" should be the same as the input Fukuoka file


class EquationSystem():
    """A class providing an interface to equation systems of all formats."""

    def __init__(self, equations, seed=0, verbose=False, order='degrevlex'):
        # initialize the instance with a list of equations over a multivariate ring
        self.seed = seed
        self.verbose = verbose
        # check that all equations have the same underlying ring
        assert len(equations) > 0
        assert len(set(eq.parent() for eq in equations)) == 1
        self.F = equations[0].parent().base_ring()
        self.q = self.F.order()
        # consider all equations over the ring with the least possible number of variables
        var_set = sorted(set().union(*[eq.variables() for eq in equations]), reverse=True)
        assert len(var_set) > 0
        self.order = order
        self.R = PolynomialRing(self.F, var_set, order=order)
        self.var_list = list(self.R.gens())
        # get rid of repeated exponents if over GF(2)
        equations = [delete_powers(eq) for eq in equations] if self.q == 2 else equations
        self.equations = [self.R(eq) for eq in equations if eq != 0]
        self.M = len(self.equations)
        self.N = len(self.var_list)
        self.ext_deg = self.F.modulus().degree()
        if is_prime(self.q):
            self.weil = None
        else:
            self.weil = self.apply_weil_descent()

    def apply_weil_descent(self):
        z = self.F.gens()[0]
        zs = [z ^ i for i in range(self.ext_deg)]
        F_base = GF(radical(self.q))
        # handle indexing better?
        weil_vars_str = ['w%s_%s' % (p1, p2) for p1, p2 in product(range(1, self.N + 1), range(self.ext_deg))]
        weil_ring = PolynomialRing(F_base, weil_vars_str, order=self.order)
        weil_vars = vector(weil_ring.gens())
        weil_parts = [weil_vars[i:i + self.ext_deg] for i in range(0, self.ext_deg * self.N, self.ext_deg)]
        weil_subs = vector([linear_combination(w, zs) for w in weil_parts])
        equations_weil = [eq_w for eq in self.equations for eq_w in weil_decomposition(eq(*weil_subs))]
        assert set(eq_w.parent() for eq_w in equations_weil) == {weil_ring}
        assert len(equations_weil) == self.M * self.ext_deg
        assert len(set().union(*[eq_w.variables() for eq_w in equations_weil])) == self.N * self.ext_deg
        return EquationSystem(equations_weil, self.seed, self.verbose, self.order)

    def to_degrevlex_str(self):
        degrevlex_mons = compute_degrevlex_mons(self.var_list)
        coeff_str = ""
        for eq in self.equations:
            coeffs = [eq.monomial_coefficient(mon) for mon in degrevlex_mons[:-1]] + [eq.constant_coefficient()]
            coeff_str += " ".join([elt_to_str(self.q, coeff) for coeff in coeffs]) + " ;\n"
        return coeff_str

    def save_one(self, eq_format, file_path):
        # handle overwrite behaviour
        if file_path.is_file() and self.verbose:
            print("The file {} already exists!".format(str(file_path)))
            return

        if eq_format == 'anf':
            var_list = self.var_list
            # assume sorted variables so that e.g. x1*x2 always appears instead of x2*x1
            var_prod_dict = {v1 * v2: sorted([i + 1, j + 1]) for i, v1 in enumerate(
                var_list) for j, v2 in enumerate(var_list) if v1 != v2}
            with open(file_path, 'w') as f:
                f.write("p anf {} {}\n".format(self.N, self.M))
                for eq in self.equations:
                    const_present = False
                    anf_line = "x "
                    for mon in eq.monomials():
                        if mon in var_list:
                            anf_line += "{} ".format(var_list.index(mon) + 1)
                        elif mon in var_prod_dict.keys():
                            anf_line += ".2 {} {} ".format(
                                var_prod_dict[mon][0], var_prod_dict[mon][1])
                        else:
                            assert mon == 1
                            const_present = True
                    if not const_present:
                        # the right hand side of the equation must correspond to True
                        anf_line += "T "
                    anf_line += "0\n"
                    f.write(anf_line)

        elif eq_format == 'cb_gpu':
            '''The format for the GPU F2 crossbred solver of Niederhagen, Ning and Yang: https://github.com/kcning/mqsolver/'''
            with open(file_path, 'w') as f:
                f.write(
                    """Galois Field : GF(2)
Number of variables (n) : {}
Number of polynomials (m) : {}
Seed : {}
Order : graded reverse lex order

*********************\n""".format(self.N, self.M, self.seed))
                f.write(self.to_degrevlex_str())

        elif eq_format == 'cb_orig':
            with open(file_path, 'w') as f:
                for eq in self.equations:
                    eq_repr = []
                    for var_tuple, coeff in eq.dict().items():
                        if coeff == 1:
                            # create an integer whose binary representation corresponds to variables present in the monomial
                            mon_repr = ZZ(''.join([str(k) for k in var_tuple]), 2)
                            eq_repr.append(mon_repr)
                    for mon_repr in sorted(eq_repr, reverse=True):
                        f.write(str(mon_repr) + "\n")
                    f.write(str(-1) + "\n")

        elif eq_format == 'cnf':
            var_list = self.var_list
            var_prod_list = []
            M = self.M
            N = self.N
            with open(file_path, 'w') as f:
                f.write("p cnf {} {}\n".format(
                    N + binomial(N, 2) + 1, M + 3 * binomial(N, 2) + 1))
                # introduce the constant variable
                f.write("{} 0\n".format(N + binomial(N, 2) + 1))
                # convert ANDs to ORs by introducing new variables
                prod_index = 0
                for i, _ in enumerate(var_list):
                    for j in range(i + 1, len(var_list)):
                        prod_index += 1
                        var_prod_list.append(var_list[i] * var_list[j])
                        f.write("{} -{} 0\n".format(i + 1, N + prod_index))
                        f.write("{} -{} 0\n".format(j + 1, N + prod_index))
                        f.write("-{} -{} {} 0\n".format(i +
                                                        1, j + 1, N + prod_index))
                for eq in self.equations:
                    const_present = False
                    cnf_line = "x "
                    for mon in eq.monomials():
                        # add linear monomials
                        if mon in var_list:
                            cnf_line += "{} ".format(var_list.index(mon) + 1)
                        # add quadratic monomials
                        elif mon in var_prod_list:
                            cnf_line += "{} ".format(
                                N + var_prod_list.index(mon) + 1)
                        else:
                            assert mon == 1
                            const_present = True
                    if not const_present:
                        # the right hand side of the equation must correspond to True
                        cnf_line += "{} ".format(str(N + binomial(N, 2) + 1))
                    cnf_line += "0\n"
                    f.write(cnf_line)

        elif eq_format == 'magma':
            q = self.q
            if q == 2:
                # for optimized performance
                ring_type = "BooleanPolynomialRing("
                field_string = f"F := GaloisField({q});\n"
            else:
                ring_type = "PolynomialRing(F, "
                field_string = f"F<{self.F.gen()}> := GaloisField({q});\n"
            var_list = [str(var) for var in self.var_list]
            with open(file_path, 'w') as f:
                f.write(field_string)
                f.write(f"R<{', '.join(var_list)}> := {ring_type}{len(var_list)}, \"grevlex\");\n")
                f.write(f"I := ideal<R |\n")
                # add field equations
                for v in var_list:
                    f.write(f"{v}^{q} - {v},\n")
                # add system equations
                for i, eq in enumerate(self.equations):
                    f.write(f"{eq}")
                    if i != len(self.equations) - 1:
                        f.write(",\n")
                    else:
                        f.write("\n")
                f.write(f">;\n")
                # use the F4 algorithm
                f.write("GroebnerBasis(I: Faugere:=true);\n")
                f.write("Variety(I);")

        elif eq_format == 'mq':
            var_list = [str(var) for var in self.var_list]
            with open(file_path, 'w') as f:
                f.write("# Variables:\n")
                f.write(', '.join(var_list) + "\n#\n")
                f.write("# Equations:\n")
                for eq in self.equations:
                    f.write(str(eq) + "\n")

        elif eq_format == 'xl':
            '''The format for the block Wiedemann XL solver of Niederhagen: http://polycephaly.org/projects/xl'''
            with open(file_path, 'w') as f:
                f.write(self.to_degrevlex_str())

        if self.verbose:
            print("Equation system written to: " + str(file_path))

    def save_all(self, folder, base_system_name):
        eq_formats = ['anf', 'cb_gpu', 'cb_orig', 'cnf', 'magma', 'mq', 'xl']
        for eq_format in eq_formats:
            # choose Weil descent for formats intended for GF(2)
            if eq_format in ['anf', 'cb_gpu', 'cb_orig', 'cnf', 'mq'] and self.weil is not None:
                eq_path = Path(folder, f"{base_system_name}_M_{self.weil.M}_N_{self.weil.N}_weil.{eq_format}")
                self.weil.save_one(eq_format, eq_path)
            else:
                eq_path = Path(folder, f"{base_system_name}_M_{self.M}_N_{self.N}.{eq_format}")
                self.save_one(eq_format, eq_path)
