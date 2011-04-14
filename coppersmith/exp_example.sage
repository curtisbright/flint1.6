M0=cop(1.5,3)
U,T=M_U(M0,10^20,90)
M1=cop(1.501,3)
T1=U*M1
print R(log(abs(M[0]),2))," should be log of modulus"
print R(log(abs(T1[0]),2))," if smaller than log of modulus things are good"

def does_U_work(U, c):
   M1=cop(c,3)
   T1 = U*M1
   print R(log(abs(M1[0]),2))," should be log of modulus"
   print R(log(abs(T1[0]),2))," if smaller than log of modulus things are good"

does_U_work(U,1.515)
