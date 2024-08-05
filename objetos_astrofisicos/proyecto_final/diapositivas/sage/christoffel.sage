m = var('M')
M = Manifold(4, 'R^4', start_index=1)

c_spher.<t,r,th,ph> = M.chart(r't:(0,+oo) r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')

g = M.metric('g')

g[1,1], g[2,2], g[3,3], g[4,4] = (-1)*(1-2*m/r), (1-2*m/r)^(-1), r^2, r^2*sin(th)^2

print(latex(g.christoffel_symbols_display(chart=c_spher)))

