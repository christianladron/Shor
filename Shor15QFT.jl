include("ShorCompuertas.jl")
using algoritmoShor
N = 15
a=7
n = int(ceil(2*log2(N)))
baseQ=big(int(big(2)^(2*n)))
#=Ψ0 = numtoQvector(xytoQstateNum(0,1,n),n)=#
#=ρ0 = Ψ0*Ψ0';=#
pflip = 0.001

#=flips = generaTransformaciones1flip(n)=#
QFT=generaQFT(n)
#=QFTinv=generaQFTinv(n)=#
#=operadorMod = generamodN(n,N,a)=#

QFT8=open("QFT8.juliaobj","w")
serialize(QFT8,QFT)
close(QFT8)
