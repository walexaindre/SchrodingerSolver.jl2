using Symbolics

@variables x,y,z,w,r,α,e,k






function calculate_component(nvar)
    F= (α/2)*(x^2+y^2+z^2+w^2)+e*(x*y+x*z+y*w+z*w)+k*(x*w+y*z)
    Fnew = substitute(F,Dict([nvar=>r]))
    Div = (Fnew-F)/(r-nvar) 
    R = simplify(Div,expand=true,threaded=true)
    simplify(R)
end

println("x: ")
display(calculate_component(x))
println("y: ")
display(calculate_component(y))
println("z: ")
display(calculate_component(z))
println("w: ")
display(calculate_component(w))