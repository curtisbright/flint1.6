import os

R=RealField(30)

def M_in(fin):
   f=file(fin,'r')
   lngstr=f.read()
   f.close()
   lngstr=lngstr.replace('\n',',')
   lngstr=lngstr.replace(' ',',')
   if (lngstr[-1]==','):
      lngstr=lngstr[:-1]
   return matrix(eval(lngstr))

def M_QR(M, prec):
   lngstr = str(M.rows())
   lngstr=lngstr.replace(' ','')
   lngstr=lngstr.replace('(','[')
   lngstr=lngstr.replace(')',']\n')
   lngstr=lngstr.replace(',',' ')
   f=file('temp','w')
   f.write(lngstr)
   f.close()
   os.system('./qrfactor '+ str(prec) +'< temp')

def M_U(M, C, prec):
   k=M.nrows()
   I=identity_matrix(k)
   ans = I.augment(C*M)
   ans = ans.LLL()
   U=ans.submatrix(0,0,k,k)
   test = U*M
   print("checking QR on U*M")
   M_QR(test, prec)
   return U,test

#gcc Householder_QR.c -o qrfactor -std=c99 -lmpir -lmpfr -lflint -Wl,-rpath,/home/andy/Desktop/andy_flint
#./qrfactor < matrix

def cop(c, alpha):
   f=file('sollya_template', 'r')
   lngstr = f.read()
   f.close()
#here we can adjust the inputs, the elementary function, the center of the interval, desirted taylor degree, target floating point precision, and pp is the error precision
   Fs='exp(x)'
   cs=str(c)
   d=7
   p=53
   pp=110

   instr='F=' + Fs + ';\nc=' + cs + ';\nd=' + str(d) + ';\np=' + str(p) + ';\npp=' + str(pp) + ';\n'

   f=file('taylor_poly.sollya', 'w')
   f.write(instr+lngstr)
   f.close()
#in future we could approximate a cap between polynomials not f and the poly
   os.system('sollya < taylor_poly.sollya > tay_pol')

   f=file('tay_pol', 'r')
   lngstr = f.read()
   f.close()

   L=lngstr.split()
   modulus = 2^int(L[0])
   xp=floor(R(log(int(L[1]))/log(2)))
   yp=int(L[2])
   P=str(d+1)+' '
   for i in range(3, len(L), 2):
      P = P + ' ' + str( sage.rings.integer.Integer( L[i] ) * 2^( int( L[i+1] ) ) )

   f=file('pol_in','w')
   f.write(P)
   f.close()

   f=file('pol_mod','w')
   f.write(str(modulus))
   f.close()
   print("modulus power = " + L[0] + " xbound, ybound " + str(xp) + ", " + str(yp))

#gcc flint_copper.c -o copper_prog -std=c99 -lmpfr -lmpir -lflint -Wl,-rpath,/home/andy/Desktop/andy_flint
#./copper_prog alpha xp yp < pol_mod
   #alpha = 3
   ypow = d*(alpha) + 1
   os.system('./copper_prog '+str(alpha)+' '+str(xp)+' '+str(yp)+' < pol_mod')
   return M_in('pre_LLL')

def M_out(M, fstr):
   lngstr = str(M.rows())
   lngstr=lngstr.replace(' ','')
   lngstr=lngstr.replace('(','[')
   lngstr=lngstr.replace(')',']\n')
   lngstr=lngstr.replace(',',' ')
   f=file(fstr,'w')
   f.write(lngstr)
   f.close()

