module TestBigEnums
using Base.Test
@enum Test111 _zerobi=BigInt(1)
@test Integer(_zerobi) == 1
end
