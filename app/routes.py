import os
from app import app
from app.models import model
from flask import render_template, request, redirect
import scipy.stats as st
import math

print(st.norm.cdf(.25))

print(st.norm.cdf(.25))

#need to take inputs for everything. 



@app.route("/list", methods=["GET", "POST"])
def list():
    if request.method=="GET":
        print("loading get")
        return render_template('form.html')
    else:
        print("loading request")
        userdata=dict(request.form)
        So=float(userdata['stockCost'])
        Xo=float(userdata['exercisePrice'])
        To=float(userdata['timeUntilExp'])
        ro=float(userdata['riskFreeInterestRate'])
        oo=float(userdata['volatility'])
        return render_template('blackScholes.html', C=model.blackScholes(
            So, Xo, To, ro, oo))

@app.route('/index')
def loadScholes():
    return render_template('/form.html')



