using Test, BSplineExtension

P1 = CDBSplinePlatform()
P2 = SparseCDBSplinePlatform()

@test AZ_A(P1,10)≈AZ_A(P2,10)
@test AZ_Z(P1,10)≈AZ_Z(P2,10)

Matrix(AZ_A(P2,10)⊗AZ_A(P2,11))

krMatrix(AZ_A(P2,10).A)

kron(AZ_A(P2,10).A,AZ_A(P2,11).A)
