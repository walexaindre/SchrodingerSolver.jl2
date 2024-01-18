using Symbolics

@variables x,y,z,w,r,α,e,k


function calculate_component(nvar)
    F= 0.5*(7*x^2+4*y^2)
    F=-F
    Fnew = substitute(F,Dict([nvar=>r]))
    Div = (Fnew-F)/(r-nvar) 
    R = simplify(Div,expand=true,threaded=true)
    simplify(R)
end

println("x: ")
display(calculate_component(x))
println("y: ")
display(calculate_component(y))