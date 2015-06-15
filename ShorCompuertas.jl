module algoritmoShor

export numtoQvector, xytoQstateNum, generaQFT, generaQFTinv, generamodN, generaTransformaciones1flip, aplicaError

function numtoQvector(num::Int,n::Int)
    sparsevec([num+1],[1],int(2^(2*n)))
end

function xytoQstateNum(x::Int,y::Int,n::Int)
    x+y*big(2)^(n)
end

function QstateNumtoxy(num::Int,n::Int)
    maxreg = 2^n -1
    y = num>>n
    x = num & (maxreg)
    (x,y)
end

function generaQFT(n::Int)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    QFT = speye(Complex128,baseQ)
    for estado in 0:(baseQ-1)
	    QFT[estado+1,estado+1]=0
        x,y = QstateNumtoxy(estado,n)
	regxencero = int(xytoQstateNum(0,y,n))
        for estadoreg1 in 0:(Q -1)
		try
		    #=numestado = xytoQstateNum(estadoreg1,y,n)=#
		    numestado = regxencero + estadoreg1
		    #=columna += exp(2*pi*im*x*estadoreg1/Q)*numtoQvector(numestado,n)/sqrt(Q)=#
		    QFT[numestado+1,estado+1] += exp(2*pi*im*x*estadoreg1/Q)/sqrt(Q)
	    catch
		    println(x,",",y,",",estadoreg1)
	    end
        end
        #QFT[1:end,estado+1]=columna
    end
           QFT
end
function generaQFTinv(n::Int)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    QFTinv = speye(Complex128,baseQ)
    for estado in 0:(baseQ-1)
	    QFTinv[estado+1,estado+1]=0
        x,y = QstateNumtoxy(estado,n)
	regxencero = int(xytoQstateNum(0,y,n))
        for estadoreg1 in 0:(Q -1)
	    numestado = regxencero + estadoreg1
	    QFTinv[numestado+1,estado+1] += exp(-2*pi*im*x*estadoreg1/Q)/sqrt(Q)
        end
    end
    QFTinv
end

function generamodN(n::Int,N::Int,a::Int)
    Q=int(big(2)^n)
    baseQ=int(Q^2)
    modN = speye(Int64,baseQ)*0
    for estado in 0:baseQ-1
	    modN[estado+1,estado+1] = 0
        x,y = QstateNumtoxy(estado,n)
        nuevoestado = y
	modestadonum = estado
            if nuevoestado < N
            nuevoestado = int((nuevoestado*big(a)^x)%N)
        modestadonum = xytoQstateNum(x,nuevoestado,n)
            end
   #     println(x,",",y,",",nuevoestado)
        modN[modestadonum+1,estado+1] = 1
    end
    modN
end

function generaTransformaciones1flip(n::Int)
    Qbits = 2*n
    dosote = big(2)
    baseQ = int(big(2)^(2*n))
    transformaciones = Array(SparseMatrixCSC{Int64,Int64},Qbits)
    for mutación in 0:int(Qbits)-1
        estado1flip = speye(Int64,baseQ)
        for estado in 0:baseQ-1
		estado1flip[estado+1,estado+1] = 0
           estadomutado = int(estado $ (dosote^mutación))
           estado1flip[estadomutado+1,estado+1] = 1
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
QFT6=open("QFT6.juliaobj")
QFT = deserialize(QFT6)
close(QFT6)

QFTinv6=open("QFTinv6.juliaobj")
QFTinv = deserialize(QFTinv6)
close(QFTinv6)


mod6a5n6=open("mod6a5n6.juliaobj")
mod6a5 = deserialize(mod6a5n6)
close(mod6a5n6)


flips6=open("flips6.juliaobj")
flips = deserialize(flips6)
close(flips6)
end
