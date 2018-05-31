#Tools for basic polynomial arithmetic and computation of intersection multiplicities
#for polynomial curves over the rationals at rational points.

#Brendan Cordy, 2018

from fractions import Fraction

class AffineCurve(object):
    #Represent a polynomial in two variables using a list of lists of coefficients.
    #The value in self.coeffs[i][j] is the coefficient on the term x^(i-j)y^(j).
    def __init__(self, coeffs):
        self.coeffs = [[Fraction(x,1) for x in y] for y in coeffs]

    @classmethod
    #Construct the monomial cx^ay^b.
    def mono(cls, c, a, b):
        mon = [[0] * k for k in range(1, a + b + 2)]
        mon[a+b][b] = c
        return AffineCurve(mon)

    def __str__(self):
        #Nicely print a term with a positive coefficient.
        def print_term(c, a, b):
            term_string = ''
            if c == 1:
                if a == 0 and b == 0:
                    term_string += '1'
                else:
                    pass
            else:
                term_string += str(c)
            if a == 0:
                pass
            elif a == 1:
                term_string += 'x'
            else:
                term_string += 'x^' + str(a)
            if b == 0:
                pass
            elif b == 1:
                term_string += 'y'
            else:
                term_string += 'y^' + str(b)
            return term_string

        poly_string = ''
        #Build the string up one term at a time.
        n = len(self.coeffs) - 1
        for (i, degree_i_coeffs) in enumerate(self.coeffs[::-1]):
            for (j, c) in enumerate(degree_i_coeffs):
                if c == 0:
                    poly_string += ''
                elif c == 1:
                    poly_string += ' + ' + print_term(1, n-i-j, j)
                elif c == -1:
                    poly_string += ' - ' + print_term(1, n-i-j, j)
                elif c < 0:
                    poly_string += ' - ' + print_term(-c, n-i-j, j)
                else:
                    poly_string += ' + ' + print_term(c, n-i-j, j)

        #Eliminate extra stuff before the first term.
        if poly_string[0:3] == ' + ':
            poly_string = poly_string[3:]
        elif poly_string[0:3] == ' - ':
            poly_string = poly_string[1:]
        return poly_string

    #Throw away unnecessary zero lists at the end of the list of coefficients.
    def trim(self):
        largest_nonzero_index = 0
        for i in range(len(self.coeffs)):
            if self.coeffs[i] != [0] * (i + 1):
                largest_nonzero_index = i
        return AffineCurve(self.coeffs[:largest_nonzero_index + 1])

    def __eq__(self, other):
        return self.trim().coeffs == other.trim().coeffs

    def degree(self):
        return len(self.trim().coeffs) - 1

    def constant_term(self):
        return self.coeffs[0][0]

    def is_identically_zero(self):
        return self.degree() == 0 and self.constant_term() == 0

    def contains(self,(x,y)):
        total = 0
        for i, deg_i_coeffs in enumerate(self.coeffs):
            for j, coeff in enumerate(deg_i_coeffs):
                #Note i - j is the power on x in the term whose coefficient
                #is self.coeffs[i][j], and j is the power of y.
                total += coeff * (x ** (i - j)) * (y ** j)
        return total == 0

    #Find the largest power of x that divides the polynomial.
    def x_multiplicity(self):
        lowest_exp = self.degree()
        for i, deg_i_coeffs in enumerate(self.coeffs):
            for j, coeff in enumerate(deg_i_coeffs):
                if coeff != 0 and i - j < lowest_exp:
                    lowest_exp = i - j
        return lowest_exp

    #Find the largest power of y that divides the polynomial.
    def y_multiplicity(self):
        lowest_exp = self.degree()
        for i, deg_i_coeffs in enumerate(self.coeffs):
            for j, coeff in enumerate(deg_i_coeffs):
                if coeff != 0 and j < lowest_exp:
                    lowest_exp = j
        return lowest_exp

    #Find the coefficient and power on x in the smallest degree term without y.
    def smallest_term_without_y(self):
        coeff = 0
        for i in range(0, self.degree() + 1):
            #Note that self.coeffs[i][0] is a degree i term without ys.
            if self.coeffs[i][0] != 0:
                return self.coeffs[i][0], i
        return 0, 0

    #Find the coefficient and power on x in the largest degree term without y.
    def largest_term_without_y(self):
        coeff = 0
        largest_x_exp = 0
        for i in range(self.degree(), -1, -1):
            #Note that self.coeffs[i][0] is a degree i term without ys.
            if self.coeffs[i][0] != 0:
                return self.coeffs[i][0], i
        return 0, 0

    #Divide the polynomial by x^n, if possible.
    def extract_xs(self, n):
        if self.x_multiplicity() < n:
            raise RuntimeError("Not divisible by x^" + str(n))
        else:
            result = [[0] * k for k in range(self.degree() + 1 - n)]
            for i in range(self.degree() + 1 - n):
                result[i] = [self.coeffs[i + n][j] for j in range(i + 1)]
            return AffineCurve(result)

    #Divide the polynomial by y^n, if possible.
    def extract_ys(self, n):
        if self.y_multiplicity() < n:
            raise RuntimeError("Not divisible by y^" + str(n))
        else:
            result = [[0] * k for k in range(self.degree() + 1 - n)]
            for i in range(self.degree() + 1 - n):
                result[i] = [self.coeffs[i + n][j + n] for j in range(i + 1)]
            return AffineCurve(result)

    def __add__(self, other):
        result_degree = max(self.degree(), other.degree())
        result = [[] for k in range(result_degree + 1)]
        for i in range(result_degree + 1):
            if i > self.degree():
                result[i] = other.coeffs[i]
            elif i > other.degree():
                result[i] = self.coeffs[i]
            else:
                result[i] = [self.coeffs[i][j] + other.coeffs[i][j] for j in range(i + 1)]
        return AffineCurve(result).trim()

    def __mul__(self, other):
        if type(other) is int or type(other) is Fraction:
            return AffineCurve([[other * x for x in y] for y in self.coeffs])
        else:
            n, m = self.degree(), other.degree()
            result_degree = n + m
            result = AffineCurve([[0] * k for k in range(1, result_degree + 1)])
            for i_1 in range(n + 1):
                for j_1 in range(i_1 + 1):
                    a = self.coeffs[i_1][j_1]
                    for i_2 in range(m + 1):
                        for j_2 in range(i_2 + 1):
                            b = other.coeffs[i_2][j_2]
                            result += AffineCurve.mono(a * b, i_1 + i_2 - j_1 - j_2, j_1 + j_2)
            return result.trim()

    def __rmul__(self, other):
        return self * other

    def pow(self, n):
        if n == 0:
            return 1
        elif n == 1:
            return self
        elif n > 1:
            result = AffineCurve(self.coeffs)
            for i in range(n - 1):
                result *= AffineCurve(self.coeffs)
            return result

    #Return the new polynomial after changing coordinates so the origin is at P.
    def shift(self, P):
        x_0, y_0 = P
        result = AffineCurve([self.coeffs[0]])
        for i in range(1, self.degree() + 1):
            for j in range(i + 1):
                x_poly = AffineCurve([[x_0],[1,0]]).pow(i - j)
                y_poly = AffineCurve([[y_0],[0,1]]).pow(j)
                shifted_term = self.coeffs[i][j] * x_poly * y_poly
                result += shifted_term
        return result

#Compute the intersection multiplicity of two affine curves at the origin.
def I_0(F, G):
    #If either curve does not pass through the origin, the multiplicity is zero.
    if F.constant_term() != 0 or G.constant_term() != 0:
        return 0
    #If one curve is identically zero, the multiplicity is infinity.
    if F.is_identically_zero() or G.is_identically_zero():
        return "Infinite"
    #If both curves pass through the origin and neither is identically zero.
    else:
        #Find the coefficient and power of the largest term of each polynomial which
        #is independent of y.
        a, m = F.largest_term_without_y()
        b, n = G.largest_term_without_y()
        #If both polynomials contain a term independent of y, say ax^m in F and bx^n
        #in G, with m >= n, I(F,G) = I(F - (b/a)x^(m-n)G, G), and F - (b/a)x^(m-n)G
        #has smaller degree than F.
        if a != 0 and b != 0:
            if m >= n:
                return I_0(F + AffineCurve.mono(Fraction(-a, b), m - n, 0) * G, G)
            else:
                return I_0(G + AffineCurve.mono(Fraction(-b, a), n - m, 0) * F, F)
        #If F doesn't have a term independent of y (F is divisible by y) but G does.
        elif a == 0 and b != 0:
            #Find the lowest powered term independent of y in G. If it's cx^k, then
            #we have I(F,G) = I(y,G) + I(F',G) = k + I(F',G), where yF' = F.
            c, k = G.smallest_term_without_y()
            return k + I_0(F.extract_ys(1), G)
        #If G doesn't have a term independent of y (G is divisible by y) but F does,
        #swap the roles of the two polynomials and use the previous clause.
        elif a != 0 and b == 0:
            return I_0(G, F)
        #If both curves don't contain a term without y, then both are divisible by y,
        #so the curves have a common component and the multiplicity is infinite.
        else:
            return "Infinite"

#Compute the intersection multiplicity of two affine curves at any point.
def I(F, G, P):
    return I_0(F.shift(P), G.shift(P))
