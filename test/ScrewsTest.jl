using Screws

m = 0.2
Iz = Matrix3(Diagonal([0.5, 0.5, 0.0375]))
r_A = Vector3(0.75,0.0,0.0)
cg = Vector3(0.45,0.0,0.0)
n = Wrench(ĵ, r_A)
v = Twist(-5.0*k̂, ô)
I = spi(m,Iz,cg)
e = 0
s = Twist(k̂, ô)
T = I*s/dot(s,I*s)
Λ = s*s'/dot(s,I*s)
J = (1+e)*dot(n,v)/dot(n,Λ*n)