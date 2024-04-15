using Oscar


# choose random parameters
b = 99
k1 = rand(-b:b)
k2 = rand(-b:b)
k3 = rand(-b:b)
k4 = rand(-b:b)
k5 = rand(-b:b)
k6 = rand(-b:b)
k7 = rand(-b:b)
k8 = rand(-b:b)
k9 = rand(-b:b)
k10 = rand(-b:b)
k11 = rand(-b:b)
k12 = rand(-b:b)
k13 = rand(-b:b)
k14 = rand(-b:b)
k15 = rand(-b:b)
k16 = rand(-b:b)
k17 = rand(-b:b)
k18 = rand(-b:b)
k19 = rand(-b:b)
k20 = rand(-b:b)
k21 = rand(-b:b)
k22 = rand(-b:b)
k23 = rand(-b:b)
k24 = rand(-b:b)
k25 = rand(-b:b)
k26 = rand(-b:b)
k27 = rand(-b:b)
k28 = rand(-b:b)
k29 = rand(-b:b)
k30 = rand(-b:b)
k31 = rand(-b:b)
c1 = rand(-b:b)
c2 = rand(-b:b)
c3 = rand(-b:b)
c4 = rand(-b:b)
c5 = rand(-b:b)

R,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17, w18, w19) = polynomial_ring(QQ,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19", "w1", "w2", "w3", "w4", "w5", "w6", "w7", "w8", "w9", "w10", "w11", "w12", "w13", "w14", "w15", "w16", "w17", "w18", "w19"])



linearSystem = [k27*w16 - (k1*k26)//k2*w18,
                k29*w8 - k28*w12 - (k6*k8)//(k7 + k8)*w13 + (k1*k3*k5)//(k2*k4 + k2*k5)*w19,
                (k17*k19)//(k18 + k19)*w9 - (k14*k16)//(k15 + k16)*w17,
                -k29*w8 - (k17*k19)//(k18 + k19)*w9 + k28*w12 + (k14*k16)//(k15 + k16)*w17,
                k12*w1 + k31*w3 + (-k13 - k30)*w5 - (k9*k11)//(k10 + k11)*w15,
                (-k23 - k31)*w3 + k30*w5 - (k20*k22)//(k21 + k22)*w11,
                -c1*w1 + w16 + k14//(k15 + k16)*w17 + (k1 + k2)//k2*w18 + (k1*k3)//(k2*k4 + k2*k5)*w19,
                -c2*w1 + w8 + k17//(k18 + k19)*w9 + w10 + k20//(k21 + k22)*w11 + w12 + k6//(k7 + k8)*w13 + w14 + k9//(k10 + k11)*w15 + k14//(k15 + k16)*w17 + (k1*k3)//(k2*k4 + k2*k5)*w19,
                -c3*w1 + w7 + k6//(k7 + k8)*w13,
                -c4*w1 + w6 + k17//(k18 + k19)*w9,
                -c5*w1 + w2 + k24//k25*w4]

conservationLaws = [
    (x1+x2+x3+x14+x15)-c1,
	  (x4+x5+x6+x7+x14+x15+x16+x17+x18+x19)-c2,
	  (x8+x16)-c3,
	  (x9+x17)-c4,
	  (x12+x13)-c5];

binomialSystem = [w1 - 1,
                  -x12 + w2,
                  -x11 + w3,
                  -x11*x12 + w4,
                  -x10 + w5,
                  -x9 + w6,
                  -x8 + w7,
                  -x7 + w8,
                  -x7*x9 + w9,
                  -x6 + w10,
                  -x6*x11 + w11,
                  -x5 + w12,
                  -x5*x8 + w13,
                  -x4 + w14,
                  -x4*x10 + w15,
                  -x3 + w16,
                  -x3*x6 + w17,
                  -x1 + w18,
                  -x1*x4 + w19]
