︠09ca892f-60a8-4394-b576-f8b4c0399f14︠
from sage.matrix.constructor import matrix

def nullspace_of_polynomialsR(S):
    # Get the monomials appearing in all polynomials
    all_monomials = set()
    for poly in S:
        all_monomials.update(poly.monomials())
    
    # Create coordinate vectors for each polynomial
    coord_vectors = []
    for polynomial in S:
        coeffs = [polynomial.monomial_coefficient(monomial) for monomial in all_monomials]
        coord_vectors.append(vector(coeffs))
    
    # Matrix with coordinate vectors as columns
    A = matrix(coord_vectors).transpose()
    
    # Print original matrix
    print("Original Matrix:")
    print(A)

    # Rank of the matrix
    rank = A.rank()
    
    print("Matrix in Reduced Row Echelon form:")
    print(A.rref())
    
    # Nullspace of the matrix
    nullspace = A.right_kernel()
    
    # Check if polynomials are linearly independent
    lin_indp = 1 if rank == len(S) else 0
    
    return lin_indp, rank, nullspace

# Example usage
R.<x,y> = PolynomialRing(ZZ, ['x', 'y'], order='degrevlex')
f1 = x^2
f2 = x^1
f3 = x^3
f4 = 4*(x^5 + y)
S = {f1, f2, f3, f4}

S1 = {f1, f2, f3}

lin_indp, rank, nullspace = nullspace_of_polynomialsR(S)
print("Rank of the matrix:", rank, "Lin indp:",lin_indp)
print("Basis for the nullspace:")
print(nullspace.basis())

lin_indp, rank, nullspace = nullspace_of_polynomialsR(S1)
print("Rank of the matrix:", rank, "Lin indp:",lin_indp)
print("Basis for the nullspace:")
print(nullspace.basis())

︡c0e015a9-5d6e-457a-a8e0-36a09e30b9a4︡{"stdout":"Original Matrix:\n[4 0 0 0]\n[4 0 0 0]\n[0 1 0 0]\n[0 0 1 0]\n[0 0 0 1]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0]\n[0 1 0 0]\n[0 0 1 0]\n[0 0 0 1]\n[0 0 0 0]\n"}︡{"stdout":"Rank of the matrix: 4 Lin indp: 1\n"}︡{"stdout":"Basis for the nullspace:\n"}︡{"stdout":"[\n\n]\n"}︡{"stdout":"Original Matrix:\n[1 0 0]\n[0 1 0]\n[0 0 1]\nMatrix in Reduced Row Echelon form:\n[1 0 0]\n[0 1 0]\n[0 0 1]\n"}︡{"stdout":"Rank of the matrix: 3 Lin indp: 1\n"}︡{"stdout":"Basis for the nullspace:\n"}︡{"stdout":"[\n\n]\n"}︡{"done":true}
︠07787935-e8dd-4086-9602-343c7b714b5bs︠
from sage.matrix.constructor import matrix

def fglm(S, G1):
    # Initialize the new Groebner basis and the set of linearly independent polynomials
    G2 = list()
    independent_polynomials = list()
    independent_polynomials.append(R.gen(1)^0)
    check_polys = independent_polynomials.copy()
    
    original_polys = list()
    original_polys.append(R.gen(1)^0)
    # Iterate through each variable
    for i in [1, 0]:

        exp = 1
        while True:
            # Reduce the monomial by G1
            original = R.gen(i)^exp
            original_polys.append(original)
            print("original polys:",original_polys)
            
            reduced = (original).reduce(G1)
            print("Val:",R.gen(i),"Processing exponent:", exp, "Reduced:", reduced)


            print("i",i,"#1",independent_polynomials)

            check_polys.append(reduced)

            print("i",i,"#1",independent_polynomials)


            print("check_polys",check_polys,"independent_polynomials",independent_polynomials)
            # Check if the reduced monomial is a linear combination of elements in independent_polynomials
            lin_indp, _, nullspace = nullspace_of_polynomialsR(check_polys)
            print("reduced:",reduced,"linindp",lin_indp)
            if lin_indp:
                independent_polynomials.append(reduced)

                exp += 1 # inc exp
            else:
                # Extract coefficients from a vector in the nullspace

                basis_matrix = nullspace.basis_matrix()
                coefficients = basis_matrix[0]  # Take the first vector in the basis
                

                print("check_polys now:",check_polys,"original_polys now:",original_polys)

                print("coefficients",basis_matrix)
                print("current polys:",independent_polynomials)
                
                dot_product = coefficients * vector(original_polys)
                G2.append(dot_product)
                
                original_polys.pop()
                check_polys.pop()

                break

    return G2

# Example input
R.<x,y> = PolynomialRing(QQ,2,order='degrevlex')
f1 = y^3+x^2
f2 = x^2*y+x^2

S = [f1, f2]

I = ideal(f1, f2)*R

# Compute the Groebner basis G1 in degrevlex order
G1 = I.groebner_basis()

# New lex order Groebner basis
G2 = fglm(S, G1)

print("G1 (degrevlex order):", G1)
print("G2 (lex order):", G2)
︡450f0f22-a4ca-469e-a996-70514260c545︡{"stdout":"original polys: [1, y]\nVal: y Processing exponent: 1 Reduced: y\ni 1 #1 [1]\ni 1 #1 [1]\ncheck_polys [1, y] independent_polynomials [1]\nOriginal Matrix:\n[1 0]\n[0 1]\nMatrix in Reduced Row Echelon form:\n[1 0]\n[0 1]\nreduced: y linindp 1\noriginal polys: [1, y, y^2]\nVal: y Processing exponent: 2 Reduced: y^2\ni 1 #1 [1, y]\ni 1 #1 [1, y]\ncheck_polys [1, y, y^2] independent_polynomials [1, y]\nOriginal Matrix:\n[1 0 0]\n[0 1 0]\n[0 0 1]\nMatrix in Reduced Row Echelon form:\n[1 0 0]\n[0 1 0]\n[0 0 1]\nreduced: y^2 linindp 1\noriginal polys: [1, y, y^2, y^3]\nVal: y Processing exponent: 3 Reduced: -x^2\ni 1 #1 [1, y, y^2]\ni 1 #1 [1, y, y^2]\ncheck_polys [1, y, y^2, -x^2] independent_polynomials [1, y, y^2]\nOriginal Matrix:\n[ 1  0  0  0]\n[ 0  1  0  0]\n[ 0  0  1  0]\n[ 0  0  0 -1]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0]\n[0 1 0 0]\n[0 0 1 0]\n[0 0 0 1]\nreduced: -x^2 linindp 1\noriginal polys: [1, y, y^2, y^3, y^4]\nVal: y Processing exponent: 4 Reduced: x^2\ni 1 #1 [1, y, y^2, -x^2]\ni 1 #1 [1, y, y^2, -x^2]\ncheck_polys [1, y, y^2, -x^2, x^2] independent_polynomials [1, y, y^2, -x^2]\nOriginal Matrix:\n[ 1  0  0  0  0]\n[ 0  1  0  0  0]\n[ 0  0  1  0  0]\n[ 0  0  0 -1  1]\nMatrix in Reduced Row Echelon form:\n[ 1  0  0  0  0]\n[ 0  1  0  0  0]\n[ 0  0  1  0  0]\n[ 0  0  0  1 -1]\nreduced: x^2 linindp 0\ncheck_polys now: [1, y, y^2, -x^2, x^2] original_polys now: [1, y, y^2, y^3, y^4]\ncoefficients [0 0 0 1 1]\ncurrent polys: [1, y, y^2, -x^2]\noriginal polys: [1, y, y^2, y^3, x]\nVal: x Processing exponent: 1 Reduced: x\ni 0 #1 [1, y, y^2, -x^2]\ni 0 #1 [1, y, y^2, -x^2]\ncheck_polys [1, y, y^2, -x^2, x] independent_polynomials [1, y, y^2, -x^2]\nOriginal Matrix:\n[ 1  0  0  0  0]\n[ 0  1  0  0  0]\n[ 0  0  1  0  0]\n[ 0  0  0 -1  0]\n[ 0  0  0  0  1]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0 0]\n[0 1 0 0 0]\n[0 0 1 0 0]\n[0 0 0 1 0]\n[0 0 0 0 1]\nreduced: x linindp 1\noriginal polys: [1, y, y^2, y^3, x, x^2]\nVal: x Processing exponent: 2 Reduced: x^2\ni 0 #1 [1, y, y^2, -x^2, x]\ni 0 #1 [1, y, y^2, -x^2, x]\ncheck_polys [1, y, y^2, -x^2, x, x^2] independent_polynomials [1, y, y^2, -x^2, x]\nOriginal Matrix:\n[ 1  0  0  0  0  0]\n[ 0  1  0  0  0  0]\n[ 0  0  1  0  0  0]\n[ 0  0  0 -1  0  1]\n[ 0  0  0  0  1  0]\nMatrix in Reduced Row Echelon form:\n[ 1  0  0  0  0  0]\n[ 0  1  0  0  0  0]\n[ 0  0  1  0  0  0]\n[ 0  0  0  1  0 -1]\n[ 0  0  0  0  1  0]\nreduced: x^2 linindp 0\ncheck_polys now: [1, y, y^2, -x^2, x, x^2] original_polys now: [1, y, y^2, y^3, x, x^2]\ncoefficients [0 0 0 1 0 1]\ncurrent polys: [1, y, y^2, -x^2, x]\n"}︡{"stdout":"G1 (degrevlex order): [x^4 - x^2, x^2*y + x^2, y^3 + x^2]\n"}︡{"stdout":"G2 (lex order): [y^4 + y^3, y^3 + x^2]\n"}︡{"done":true}
︠497fa1a2-b72f-41d9-b290-d1525f04ecc5︠
from sage.matrix.constructor import matrix

def fglm(S, G1):
    # Initialize the new Groebner basis and the set of linearly independent polynomials
    G2 = list()
    independent_polynomials = list()
    independent_polynomials.append(R.gen(1)^0)
    check_polys = independent_polynomials.copy()
    
    original_polys = list()
    original_polys.append(R.gen(1)^0)
    # Iterate through each variable

    for i in [2, 1, 0]:  

        exp = 1
        while True:
            # Reduce the monomial by G1
            original = R.gen(i)^exp
            original_polys.append(original)
            print("original polys:",original_polys)
            
            reduced = (original).reduce(G1)
            print("Val:",R.gen(i),"Processing exponent:", exp, "Reduced:", reduced)

            print("i",i,"#1",independent_polynomials)

            check_polys.append(reduced)

            print("i",i,"#1",independent_polynomials)


            print("check_polys",check_polys,"independent_polynomials",independent_polynomials)
            # Check if the reduced monomial is a linear combination of elements in independent_polynomials
            lin_indp, _, nullspace = nullspace_of_polynomialsR(check_polys)
            print("reduced:",reduced,"linindp",lin_indp)
            if lin_indp:
                independent_polynomials.append(reduced)

                exp += 1 # inc exp
            else:
                # Extract coefficients from a vector in the nullspace

                basis_matrix = nullspace.basis_matrix()
                coefficients = basis_matrix[0]  # Take the first vector in the basis
                

                print("check_polys now:",check_polys,"original_polys now:",original_polys)

                print("coefficients",basis_matrix)
                print("current polys:",independent_polynomials)
                
                dot_product = coefficients * vector(original_polys)
                G2.append(dot_product)
                
                print("G2 is", G2)
                
                original_polys.pop()
                check_polys.pop()

                break
 
    return G2

R.<z,y,x> = PolynomialRing(QQ,3,order='degrevlex')
f1 = y^3 
f2 = x^2*y + x^2 + y + y^3
f3 = z*5 + y^3

S = [f1, f2, f3]

I = ideal(f1, f2, f3)*R

# Compute the Groebner basis G1 in degrevlex order
G1 = I.groebner_basis()
print("Original Groebner basis", G1)

# New lex order Groebner basis
G2 = fglm(S, G1)

print("G1 (degrevlex order):", G1)
print("G2 (lex order):", G2)
︡e7691148-9dd3-433d-924d-1cd88b123f20︡{"stdout":"Original Groebner basis [x^4 - x^2 - y, y*x^2 + x^2 + y, y^2 - x^2 - y, z]\n"}︡{"stdout":"original polys: [1, x]\nVal: x Processing exponent: 1 Reduced: x\ni 2 #1 [1]\ni 2 #1 [1]\ncheck_polys [1, x] independent_polynomials [1]\nOriginal Matrix:\n[1 0]\n[0 1]\nMatrix in Reduced Row Echelon form:\n[1 0]\n[0 1]\nreduced: x linindp 1\noriginal polys: [1, x, x^2]\nVal: x Processing exponent: 2 Reduced: x^2\ni 2 #1 [1, x]\ni 2 #1 [1, x]\ncheck_polys [1, x, x^2] independent_polynomials [1, x]\nOriginal Matrix:\n[1 0 0]\n[0 0 1]\n[0 1 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0]\n[0 1 0]\n[0 0 1]\nreduced: x^2 linindp 1\noriginal polys: [1, x, x^2, x^3]\nVal: x Processing exponent: 3 Reduced: x^3\ni 2 #1 [1, x, x^2]\ni 2 #1 [1, x, x^2]\ncheck_polys [1, x, x^2, x^3] independent_polynomials [1, x, x^2]\nOriginal Matrix:\n[1 0 0 0]\n[0 0 1 0]\n[0 0 0 1]\n[0 1 0 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0]\n[0 1 0 0]\n[0 0 1 0]\n[0 0 0 1]\nreduced: x^3 linindp 1\noriginal polys: [1, x, x^2, x^3, x^4]\nVal: x Processing exponent: 4 Reduced: x^2 + y\ni 2 #1 [1, x, x^2, x^3]\ni 2 #1 [1, x, x^2, x^3]\ncheck_polys [1, x, x^2, x^3, x^2 + y] independent_polynomials [1, x, x^2, x^3]\nOriginal Matrix:\n[1 0 0 0 0]\n[0 0 0 0 1]\n[0 0 1 0 1]\n[0 0 0 1 0]\n[0 1 0 0 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0 0]\n[0 1 0 0 0]\n[0 0 1 0 0]\n[0 0 0 1 0]\n[0 0 0 0 1]\nreduced: x^2 + y linindp 1\noriginal polys: [1, x, x^2, x^3, x^4, x^5]\nVal: x Processing exponent: 5 Reduced: x^3 + y*x\ni 2 #1 [1, x, x^2, x^3, x^2 + y]\ni 2 #1 [1, x, x^2, x^3, x^2 + y]\ncheck_polys [1, x, x^2, x^3, x^2 + y, x^3 + y*x] independent_polynomials [1, x, x^2, x^3, x^2 + y]\nOriginal Matrix:\n[1 0 0 0 0 0]\n[0 0 0 0 1 0]\n[0 0 0 0 0 1]\n[0 0 1 0 1 0]\n[0 0 0 1 0 1]\n[0 1 0 0 0 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0 0 0]\n[0 1 0 0 0 0]\n[0 0 1 0 0 0]\n[0 0 0 1 0 0]\n[0 0 0 0 1 0]\n[0 0 0 0 0 1]\nreduced: x^3 + y*x linindp 1\noriginal polys: [1, x, x^2, x^3, x^4, x^5, x^6]\nVal: x Processing exponent: 6 Reduced: 0\ni 2 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ni 2 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ncheck_polys [1, x, x^2, x^3, x^2 + y, x^3 + y*x, 0] independent_polynomials [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nOriginal Matrix:\n[1 0 0 0 0 0 0]\n[0 0 0 0 1 0 0]\n[0 0 0 0 0 1 0]\n[0 0 1 0 1 0 0]\n[0 0 0 1 0 1 0]\n[0 1 0 0 0 0 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0 0 0 0]\n[0 1 0 0 0 0 0]\n[0 0 1 0 0 0 0]\n[0 0 0 1 0 0 0]\n[0 0 0 0 1 0 0]\n[0 0 0 0 0 1 0]\nreduced: 0 linindp 0\ncheck_polys now: [1, x, x^2, x^3, x^2 + y, x^3 + y*x, 0] original_polys now: [1, x, x^2, x^3, x^4, x^5, x^6]\ncoefficients [0 0 0 0 0 0 1]\ncurrent polys: [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nG2 is [x^6]\noriginal polys: [1, x, x^2, x^3, x^4, x^5, y]\nVal: y Processing exponent: 1 Reduced: y\ni 1 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ni 1 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ncheck_polys [1, x, x^2, x^3, x^2 + y, x^3 + y*x, y] independent_polynomials [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nOriginal Matrix:\n[1 0 0 0 0 0 0]\n[0 0 0 0 1 0 1]\n[0 0 0 0 0 1 0]\n[0 0 1 0 1 0 0]\n[0 0 0 1 0 1 0]\n[0 1 0 0 0 0 0]\nMatrix in Reduced Row Echelon form:\n[ 1  0  0  0  0  0  0]\n[ 0  1  0  0  0  0  0]\n[ 0  0  1  0  0  0 -1]\n[ 0  0  0  1  0  0  0]\n[ 0  0  0  0  1  0  1]\n[ 0  0  0  0  0  1  0]\nreduced: y linindp 0\ncheck_polys now: [1, x, x^2, x^3, x^2 + y, x^3 + y*x, y] original_polys now: [1, x, x^2, x^3, x^4, x^5, y]\ncoefficients [ 0  0  1  0 -1  0  1]\ncurrent polys: [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nG2 is [x^6, -x^4 + x^2 + y]\noriginal polys: [1, x, x^2, x^3, x^4, x^5, z]\nVal: z Processing exponent: 1 Reduced: 0\ni 0 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ni 0 #1 [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\ncheck_polys [1, x, x^2, x^3, x^2 + y, x^3 + y*x, 0] independent_polynomials [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nOriginal Matrix:\n[1 0 0 0 0 0 0]\n[0 0 0 0 1 0 0]\n[0 0 0 0 0 1 0]\n[0 0 1 0 1 0 0]\n[0 0 0 1 0 1 0]\n[0 1 0 0 0 0 0]\nMatrix in Reduced Row Echelon form:\n[1 0 0 0 0 0 0]\n[0 1 0 0 0 0 0]\n[0 0 1 0 0 0 0]\n[0 0 0 1 0 0 0]\n[0 0 0 0 1 0 0]\n[0 0 0 0 0 1 0]\nreduced: 0 linindp 0\ncheck_polys now: [1, x, x^2, x^3, x^2 + y, x^3 + y*x, 0] original_polys now: [1, x, x^2, x^3, x^4, x^5, z]\ncoefficients [0 0 0 0 0 0 1]\ncurrent polys: [1, x, x^2, x^3, x^2 + y, x^3 + y*x]\nG2 is [x^6, -x^4 + x^2 + y, z]"}︡{"stdout":"\n"}︡{"stdout":"G1 (degrevlex order): [x^4 - x^2 - y, y*x^2 + x^2 + y, y^2 - x^2 - y, z]\n"}︡{"stdout":"G2 (lex order): [x^6, -x^4 + x^2 + y, z]\n"}︡{"done":true}
︠c651b905-dea6-4c87-9eb9-5ef93442f4d3︠









