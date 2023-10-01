using Symbolics

@variables x,y,z,w,r,α,e,k
F= (α/2)*(x^2+y^2+z^2+w^2)+e*(x*y+x*z+y*w+z*w)+k*(x*w+y*z)

Fnew = substitute(F,Dict([w=>r]))

Div = (Fnew-F)/(r-w) 

R = simplify(Div,expand=true,threaded=true)
simplify(R)