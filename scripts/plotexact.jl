using GLMakie

xv = collect(0:0.05:2*pi);
yv = xv;
t = collect(0:0.1:10);
f(x,y,t)=exp(1im*(2x+y-t))
g(x,y,t)=exp(1im*(2x+2y-t))

#surface(,
#    colormap = :darkterrain,
#    colorrange = (80, 190),
#    axis=(type=Axis3, azimuth = pi/4))

xpoints = repeat(xv,inner=size(yv))
ypoints = repeat(yv,outer=size(xv))

zpoints = abs2.(f.(xpoints,ypoints,0))+abs2.(g.(xpoints,ypoints,0))
zpoints = Observable(zpoints)

function anim_step(t)
    zpoints[] = abs2.(f.(xpoints,ypoints,t))+abs2.(g.(xpoints,ypoints,0))
    
end

r  = surface(xpoints,ypoints,zpoints,
    colormap = :darkterrain,
    colorrange = (80, 190),
    axis=(type=Axis3, azimuth = pi/4))
display(r)
for i in t

    anim_step(i)
    r.axis.title="Step $i"
    sleep(0.1)
end