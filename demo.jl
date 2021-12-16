
# glcPointsLocal(N::Int) = [0.5(1-cos((i-1)pi/(N-1))) for i in 1:N]
# glcPointsGlobal(left::Number,right::Number,localPnts::Vector) = (right-left).*localPnts.-left

function glcPointsGlobal(left::Number,right::Number,N::Int)
    glcPointsLocal= [0.5(1-cos((i-1)pi/(N-1))) for i in 1:N]
    (right-left).*glcPointsLocal.-left
end

function numstring(x::Vector)::Float64
    pop!(cumprod(x))
end

function CofDQM(x::Vector,order::Int)
    #x: vector of sampling points;
    #order: largest derivative order;
    #return: matrix with dimension of (x,x,m)

    n = length(x)
    C = fill(0.0,(n,n,order))
    (n>order) || error("need more sampling points")

    # 1st order
    
    for i in 1:n, j in 1:n
        if i != j
            C[i,j,1] = numstring([((m!=i && m!=j) ? x[i]-x[m] : 1) for  m in 1:n])/numstring([((m!=j) ? x[j]-x[m] : 1) for  m in 1:n])

            # m_range_denominater = filter!(xx->xx!=i,[m for m in 1:n])
            # m_range_numerator = filter!(xx->xx!=j,m_range_denominater)
            #C[i,j,1] = numstring(filter(xx->xx!=0,[x[i]-x[m] for m in m_range_numerator]))/numstring(filter(xx->xx!=0,[x[j]-x[m] for m in m_range_denominater]))
            # C[i,j,1] = pop!(cumprod(filter(xx->xx!=0,[x[i]-x[m] for m in m_range_numerator])))/pop!(cumprod(filter(xx->xx!=0,[x[j]-x[m] for m in 1:n])))
        else

            # C[i,i,1] = sum(filter(xx->xx!=Inf,[1/(x[i]-x[m]) for m in 1:n]))
            C[i,i,1] = sum([((m!=i) ? 1/(x[i]-x[m]) : 0) for m in 1:n]) #这个慢了4个微秒

        
        end
    end


    # higher order
    for i in 1:n, j in 1:n, k in 2:order
        if i != j
            C[i,j,k] = k*(C[i,i,k-1]*C[i,j,1] - C[i,j,k-1]/(x[i]-x[j]))
        end
    end
    for i in 1:n, j in 1:n, k in 2:order
        if i == j
            C[i,j,k] = -sum([((m!=i) ? C[i,m,k] : 0) for m in 1:n])
        end
    end
    C
end

`
demoooo: test 1st-order approximation accuracy
`
pnts = glcPointsGlobal(-10,1,100)
# pnts = collect(0:2:10)

Coef = CofDQM(pnts,2)
f4(x::Number) = 4x^3+3x-2
f4d(x) = 12x^2+3

f4d.(pnts)-Coef[:,:,1]*f4.(pnts)

