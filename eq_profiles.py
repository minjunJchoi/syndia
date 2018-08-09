# otbain 2D profiles
import geqdsk_dk
geqdsk_fn = "../data/g013728.003900"
g1 = geqdsk_dk(filename=geqdsk_fn)
R = 2.0
z = 0.3

B = g1.B_abs(R,z)
print(B)
