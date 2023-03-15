round2 = function(x, n) {
        posneg = sign(x)
        z = abs(x)*10^n
        z = z + 0.5 + sqrt(.Machine$double.eps)
        z = trunc(z)
        z = z/10^n
        z*posneg
}
