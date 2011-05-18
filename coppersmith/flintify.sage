import os

R=RealField(30)

#*****************grand explanation comments ************************
#to make things work we need to compile the two C programs as in the comments below
#then using sage from the directory with the c programs and sollya_template (the coppersmith folder)
#run sage
#from within sage do: load 'flintify.sage'

#1) to change the function alter the line:

#Fs='exp(x)'

#and make a string version of the function you are interested in with x as the variable (the program will normalize it for us)

#2) to change the degree of the approximation change the line:

#d=7

#and make it d=5 or d=20 or whichever degree

#3) to change the precisions needed it is the lines p=53, pp=110 to whichever combination you would like

#4) To get a coppersmith matrix for the function at center point c, and checking with alpha value, alpha, 

#we use the sage command:

#M = cop( c, alpha)

#for instance M=cop(1.2, 3)

#this will actually create the coppersmith lattice (in the file pre_LLL) and reduce it (to the file #post_LLL) and load into sage the pre_LLL matrix as the output

#it also displays powers of two which give the modulus, the bound on X and the bound on Y in the #coppersmith lattice

#So now M is a coppersmith lattice

#To find a U, unimodular which reduces M the coppersmith lattice we use

#M_U(M, a constant, a precision)

#the output will be U then U*M (the reduced lattice)

#so for instance the command U, T = M_U(M, 10^20, 90) will give you the reducing U, the reduced M as T and it will display some reducedness data from T (that's what the precision does).

#We could compare the size of the first row of T and first row of M

#In SAGE to make a floating point approximation with n bits I use the following
#R=RealField(n)
#then R( log( abs( M[0]), 2 ) ) will display the bit-length of the first row of M before LLL
#and R( log( abs( T[0]), 2 ) ) the bit length of the first row of T

#Now to test if the same U will work on a nearby lattice we could make a second lattice:
#M2 = cop(1.201,3) (be sure to use the same alpha)
#test if U works on M2 by something like:

#T2 = U*M2

#then

#R( log( abs( T2[0]), 2 ) ) the bit length of the first row of T2
#the command M_U uses fpLLL from within SAGE which sometimes hangs if the dimension is too large, in reality that doesn't mean that the lattice is too big as my LLL performs better in such cases. 

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

def copx(c, alpha, xb):
   f=file('sollya_template', 'r')
   lngstr = f.read()
   f.close()
#here we can adjust the inputs, the elementary function, the center of the interval, desirted taylor degree, target floating point precision, and pp is the error precision
   Fs='exp(x)'
   cs=str(c)
   d=11
   p=53
   pp=90
#so 2^p is the target floating point precision, and pp is the range of value we are checking for in coppersmith's method
# we need to find a polynomial such that infnorm of f - P in a small interval around c will be smaller than 
# 2^pp, then we consider all intgers/2^p.  24 and 40 would make a small-scale version of TaMaDi, 53,110 would be the double precision table makers problem 
   instr='F=' + Fs + ';\nc=' + cs + ';\nd=' + str(d) + ';\np=' + str(p) + ';\npp=' + str(pp) + ';\n'

   f=file('taylor_poly.sollya', 'w')
   f.write(instr+lngstr)
   f.close()
#in future we could approximate a cap between polynomials not f and the poly
   os.system('sollya < taylor_poly.sollya > output')

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
   print("modulus power = " + L[0] + " xbound, ybound " + str(xb) + ", " + str(yp))

#gcc flint_copper.c -o copper_prog -std=c99 -lmpfr -lmpir -lflint -Wl,-rpath,/home/andy/Desktop/andy_flint
#./copper_prog alpha xp yp < pol_mod
   #alpha = 3
   print "ypow is",d*(alpha) + 1
   os.system('./copper_prog '+str(alpha)+' '+str(xb)+' '+str(yp)+' < pol_mod')
   return M_in('pre_LLL')

#copx_auto returns all the pertinent data
#M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,xb, d)
def copx_auto(c, alpha, xb, d, p, pp):
   f=file('sollya_template', 'r')
   lngstr = f.read()
   f.close()
#here we can adjust the inputs, the elementary function, the center of the interval, desirted taylor degree, target floating point precision, and pp is the error precision
   Fs='exp(x)'
   cs=str(c)
#so 2^p is the target floating point precision, and pp is the range of value we are checking for in coppersmith's method
# we need to find a polynomial such that infnorm of f - P in a small interval around c will be smaller than 
# 2^pp, then we consider all intgers/2^p.  24 and 40 would make a small-scale version of TaMaDi, 53,110 would be the double precision table makers problem 
   instr='F=' + Fs + ';\nc=' + cs + ';\nd=' + str(d) + ';\np=' + str(p) + ';\npp=' + str(pp) + ';\nxb='+str(xb)+';\n'

   f=file('taylor_poly.sollya', 'w')
   f.write(instr+lngstr)
   f.close()
#in future we could approximate a cap between polynomials not f and the poly
   os.system('sollya < taylor_poly.sollya > output')

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
   print("modulus power = " + L[0] + " xbound, ybound " + str(xb) + ", " + str(yp))

#gcc flint_copper.c -o copper_prog -std=c99 -lmpfr -lmpir -lflint -Wl,-rpath,/home/andy/Desktop/andy_flint
#./copper_prog alpha xp yp < pol_mod
   #alpha = 3
   print "ypow is",d*(alpha) + 1
   os.system('./copper_prog '+str(alpha)+' '+str(xb)+' '+str(yp)+' < pol_mod')
   return M_in('pre_LLL'), 2^int(L[0]), alpha, 2^xb, 2^yp, d*alpha + 1
#M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,xb,d)

def M_out(M, fstr):
   lngstr = str(M.rows())
   lngstr=lngstr.replace(' ','')
   lngstr=lngstr.replace('(','[')
   lngstr=lngstr.replace(')',']\n')
   lngstr=lngstr.replace(',',' ')
   f=file(fstr,'w')
   f.write(lngstr)
   f.close()

S.<x,y>=PolynomialRing(ZZ)
def poly_in(fin, ypow):
   f=file(fin,'r')
   lngstr=f.read()
   f.close()
   k = lngstr.find('  ')
   lngstr=lngstr[k+2:]
   L=list(eval(lngstr.replace(' ',',')))
   P=0
   for i in range(len(L)):
      P = P + L[i]*x^(i%ypow)*y^(floor(i/ypow))
   return P

def poly_row(M, i, ypow):
   P=0
   for j in range(M.ncols()):
      P = P + M[i][j]*x^(j%ypow)*y^(floor(j/ypow))
   return P

R=RealField(200)
def test_M(M,ypow,X,Y):
   T=M.LLL()
   print R(abs(T[0])), "should be smaller than:\n", R(abs(M[0]))
   p1 = poly_row(T,0,ypow)
   p1 = S(p1(x/X,y/Y))
   p2 = poly_row(T,1,ypow)
   p2 = S(p2(x/X,y/Y))
   f=p1.resultant(p2,y)
   f=f.univariate_polynomial()
   return f.factor()

def test_T(T,ypow,X,Y,i,j):
   print R(abs(T[0])), "should be smaller than:\n", R(abs(M[0]))
   p1 = poly_row(T,i,ypow)
   p1 = S(p1(x/X,y/Y))
   p2 = poly_row(T,j,ypow)
   p2 = S(p2(x/X,y/Y))
   f=p1.resultant(p2,y)
   f=f.univariate_polynomial()
   return f.factor()

def rec_pp_finder(alpha,d,pp,delta):
   #assume that works for pp and not pp + delta * 2
   failed = 0
   try:
      M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,5,d,pp + delta)
      pp = pp + delta
   except:
      failed = 1
   if (delta == 1):
      return pp
   else:
      return rec_pp_finder(alpha, d, pp, delta/2)

def pp_finder(c,alpha,d):
   pp = 53
   delta = 256
   failed = 0
   try:
      M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,5,d,pp)
   except:
      failed = 1
   if (failed == 1):
      return -1
   try:
      M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,5,d,pp + delta)
      failed = 2
   except:
      failed = 1
   if (failed == 2):
      return (pp + delta)
   return rec_pp_finder(alpha, d, pp, delta/2)

def lat_test(c,alpha,xb,d,pp):
   M,modulus,alpha,X,Y,ypow=copx_auto(c,alpha,xb,d,pp)
   T=M_in('post_LLL');
   return (R(abs(M[0])/abs(T[0])) > 1.01)

def rec_max_interval(c,alpha,d,pp,xb,delta):
   #assume that xb works xb + 2*delta fails
   print 'testing',xb + delta
   failed = 0
   if lat_test(c,alpha,xb + delta,d,pp):
      xb = xb + delta
   else:
      failed = 1
   if (delta == 1):
      return xb
   else:
      return rec_max_interval(c,alpha,d,pp,xb,delta/2)

def max_interval(c,alpha,d,pp):
   xb = 0
   delta = 64
   if not lat_test(c,alpha,xb,d,pp):
      return -1
   if lat_test(c,alpha,xb + delta, d, pp):
      return xb + delta
   return rec_max_interval(c,alpha,d,pp,xb,delta/2)

def alpha_d_test(c, alpha, d):
   pp = pp_finder(c,alpha,d)
   if (pp < 0):
      return -1,-1
   xb = max_interval(c,alpha,d,pp)
   return pp,xb

def msb(M, k):
   r=M.nrows()
   c=M.ncols()
   I=identity_matrix(r)
   L=[]
   mden = ceil(R(log(M.height(),2)))
   for i in range(r):
      V=M.submatrix(i,0,1,c)
      I[i,i]=ceil(R(log(V.height(),2)))-k
      if (I[i,i] < mden):
         mden = I[i,i]
      L.extend([floor(j/2^(I[i,i])) for j in V.list()])
   J=identity_matrix(r)
   for i in range(r):
      J[i,i]=2^(I[i,i]-mden)
   msbM=matrix(r,c,L)
   return I,msbM,J,mden

#the three intuition building lines at the moment
#M,modulus,alpha,X,Y,ypow=copx_auto((2^53 + 2^32)/2^53,4,32,12,300)
#rI,rowM1,rJ,mden=msb(M,250);P=rJ*rowM1;I=identity_matrix(P.nrows());IP=I.augment(P)
#TP=IP.LLL(delta=.75,eta=.81,algorithm='fpLLL:wrapper');U=TP.submatrix(0,0,TP.nrows(),TP.nrows());T=U*M;R(abs(T[0]));R(abs(M[0]))

