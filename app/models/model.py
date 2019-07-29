import scipy.stats as st
import math

def blackScholes(S, X, T, r, o):
    d1=(math.log(S/X)+(T*(r+(o**2)/2)))/(o*(T**(.5)))
    d2=d1-(o*(T**(.5)))
    C=(S*(st.norm.cdf(d1)))-(X*(math.exp(-r*T))*(st.norm.cdf(d2)))
    return C