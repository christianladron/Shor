module algoritmoShor

export numtoQvector, xytoQstateNum, generaQFT, generaQFTinv, generamodN, generaTransformaciones1flip, aplicaError

function numtoQvector(num,n)
    sparse(vcat(zeros(int(num)),[1],zeros(int(big(2)^(2*n)-num-1))))
end

function xytoQstateNum(x,y,n)
    x+y*big(2)^(n)
end

function QstateNumtoxy(num,n)
    maxreg = 2^n -1
    y = num>>n
    x = num & (maxreg)
    (x,y)
end

function generaQFT(n)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    QFT = zeros(Complex128,(baseQ,baseQ))
    for estado in 0:(baseQ-1)
        x,y = QstateNumtoxy(estado,n)
        columna = zeros(baseQ)
        for estadoreg1 in 0:(Q -1)
            numestado = xytoQstateNum(estadoreg1,y,n)
            columna += exp(2*pi*im*x*estadoreg1/Q)*numtoQvector(numestado,n)/sqrt(Q)
        end
        QFT[1:end,estado+1]=columna
    end
           sparse(QFT)
end
function generaQFTinv(n)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    QFTinv = zeros(Complex128,(baseQ,baseQ))
    for estado in 0:(baseQ-1)
        x,y = QstateNumtoxy(estado,n)
        columna = zeros(baseQ)
        for estadoreg1 in 0:(Q -1)
            numestado = xytoQstateNum(estadoreg1,y,n)
            columna += exp(-2*pi*im*x*estadoreg1/Q)*numtoQvector(numestado,n)/sqrt(Q)
        end
        QFTinv[1:end,estado+1]=columna
    end
    sparse(QFTinv)
end

function generamodN(n,N,a)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    modN = speye(baseQ,baseQ)*0
    for estado in 0:baseQ-1
        x,y = QstateNumtoxy(estado,n)
        nuevoestado = y
            if nuevoestado < N
            nuevoestado = (nuevoestado*big(a)^x)%N
            end
   #     println(x,",",y,",",nuevoestado)
        modestadonum = xytoQstateNum(x,nuevoestado,n)
        columna = numtoQvector(modestadonum,n)
        modN[:,estado+1] = columna
    end
    modN
end

function generaTransformaciones1flip(n)
    Qbits = 2*n
    dosote = big(2)
    baseQ = int(big(2)^(2*n))
    transformaciones = Array(SparseMatrixCSC{Float64,Int64},Qbits)
    for mutación in 0:int(Qbits)-1
        estado1flip = speye(baseQ)
        for estado in 0:baseQ-1
           estadomutado = estado $ (dosote^mutación)
           mutadovector = numtoQvector(estadomutado,n)
           columna = mutadovector
           estado1flip[:,estado+1] = columna
        end
        transformaciones[mutación+1] = estado1flip
    end
    transformaciones
end

function aplicaError(ρ,pmuta,transformaciones)
    dosote = big(2)
    pnomuta = 1-pmuta
    baseQ=size(ρ)[1]
    Q = int(sqrt(baseQ))
    n = int(log2(Q))
    pmuta = pmuta/(2*n)

    ρf = pnomuta*ρ
    for mutación in 1:2*n
        ρf+=pmuta*transformaciones[mutación]*ρ*transformaciones[mutación]'
    end
    ρf
end

end
