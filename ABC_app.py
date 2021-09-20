import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import dash.dependencies as dependencies
import os
import json
import plotly
import rpy2.robjects as robjects
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# from mca_lib import get_simulation_samples,counts_
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# with open('./custom.geo.json', 'r') as f:
#     data = json.load(f)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

env_num=500
import numpy as np

robjects.r(
    '''
library(magrittr)
library(dplyr)
library(spatstat)

abc_simu <-function(ndose,cidx,tys,tns,phi){
    # p.true <- c(0.05, 0.06, 0.08, 0.11, 0.19, 0.34)
    add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01)
    # phi <-0.2
    gen.u.rand <- function(k, K=5, phi=0.3, delta=0.1){
        #cps <- c(phi)
        if (k==(K+1)){
                cps <- runif(K, 0, max=phi-1*delta)
        }else if (k > 0) {
            cps <- c(runif(1, phi-1*delta, phi+1*delta))
            #cps <- c(phi)
            if (k > 1){
                cps <- c(cps, runif(k-1, min=0, max=phi-1*delta))
            }
            if (k < K){
                cps <- c(cps, runif(K-k, min=phi+1*delta, max=2*phi))
            }
        }else if (k==0){
                cps <- runif(K, min=phi+1*delta, max=2*phi)
        }
        sort(cps)
    }
    
    # generate the scenarios for Pr(Yn|M_n, Ak)
    gen.mu.rand <- function(k, J, K=5, phi=0.3, delta=0.1){
        pss <- lapply(1:J, function(i)gen.u.rand(k, K, phi, delta))
        pssMat <- do.call(rbind, pss)
        pssMat
    }
    
    gen.prior <- function(K, phi, J=1e3, delta=0.05){
        pss <- lapply(0:K, function(k)gen.mu.rand(k, J=J*(1+as.numeric(k==-1)), K=K, phi=phi, delta=delta))
        pss.prior <- do.call(rbind, pss)
        #pss.prior <- t(apply(matrix(runif(K*J, 0, 2*phi), ncol=K), 1, sort))
        pss.prior
    }
    
    earlystop <- 0
    # ndose <- length(p.true)
    # cidx <- 1#init.level
    # tys <- rep(0, ndose) # number of responses for different doses.
    # tns <- rep(0, ndose) # number of subject for different doses.
    # tys[1]<-1
    # tns[1]<-1
    
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    ps.name <- paste0("./pssprior-ndose-", ndose, "-phi-", 100*phi, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
        if (file.exists(ps.name)){
            load(ps.name)
        }else{
            pss.prior <- gen.prior(ndose, phi=phi, J=add.args$J, delta=add.args$delta)
            save(pss.prior, file=ps.name)
        }
    
    kpws.fn <- function(pss.prior, tys, tns, h=0.01){
        K <- length(tys)
        Num <- dim(pss.prior)[1]
        pss.prior.vec <- as.vector(t(pss.prior))
        tns.vec <- rep(tns, Num)
        tys.gen <- rbinom(length(tns.vec), tns.vec, pss.prior.vec)
        tys.vec <- rep(tys, Num)
        tns.mat <- matrix(tns.vec, ncol=K, byrow=T)
        diff.mat <- matrix(tys.vec - tys.gen, ncol=K, byrow=T)
        tns.mat[tns.mat==0] <- 0.1
        rate.diff.mat <- diff.mat / tns.mat
        if (is.null(h))
            h <- 0.01
        ws <- exp(-rowSums(rate.diff.mat**2)/h)# kernel fn exp(-x^2/h), rigirously, it should be exp(-x^2/2/h^2)
        ws 
    }
    
    cMTD.fn <- function(phi, pss.prior, kp.ws){
               ndose <- dim(pss.prior)[2]
               post.ms <- sapply(1:ndose, function(i)weighted.median(pss.prior[, i], w=kp.ws))
               cMTD <- which.min(abs(post.ms-phi))
               cMTD
            }
    # posterior probability of pj >= phi given data
    post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
        alp <- alp.prior + y 
        bet <- bet.prior + n - y
        1 - pbeta(phi, alp, bet)
    }
    
    overdose.fn <- function(phi, add.args=list()){
    cutoff.eli <- add.args$cutoff.eli
    cutoff.num <- add.args$cutoff.num
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    if ((pp >= cutoff.eli) & (n>=cutoff.num)){
        return(TRUE)
    }else{
        return(FALSE)
    }
    }
    cy <- tys[cidx]
    cn <- tns[cidx]
    add.args <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
    if (overdose.fn(phi, add.args)){
       tover.doses[cidx:ndose] <- 1
    }
    if (tover.doses[1] == 1){
      earlystop <- 1
       # break()
      cMTD=99
      return(cMTD)
    }
    kp.ws <- kpws.fn(pss.prior, tys, tns, h=add.args$h)
    cMTD <- cMTD.fn(phi, pss.prior, kp.ws)
    return(cMTD)
}
    '''
)

app.layout = html.Div(children=[
    html.H1(children='Approximate Bayesian Computation Design for Phase I Clinical Trials'),
    html.Div(
        children=[
            html.Label('Target:',style={'display': 'inline-block'}),
            dcc.Input(id='target', value='0.3',type='text',style={'display': 'inline-block'}),
            html.Label('Number of dose levels:',style={'display': 'inline-block'}),
            dcc.Input(id='n_doses',value=5,type='number',min=2,style={'display': 'inline-block'}),
        ],),

    html.Div([
        html.Button("Add a cohort of patients", id="add-filter", n_clicks=0),
        html.Label('Selected dose level for each cohort:',),
        html.Div(id='input-dose-container', children=[],),
        html.Label('Assigned patients number for each Cohort:',),
        html.Div(id='input-patients-container', children=[]),
        html.Label('Observed DLT number from each cohort:',),
        html.Div(id='input-dlt-container', children=[]),
        # html.Button(children='Show',id='show',n_clicks=0),
        html.Button(children='Run Simulation',id='mc',n_clicks=0),
        html.Div([
                    html.Label(children='Current MTD: ',style={'display': 'inline-block'}),
                    html.Output(id='cMTD',children='0',style={'display': 'inline-block'}),
                ]),
            html.Div([
                        html.Label(children='Suggestion for Dose Level of Next Cohort: ',style={'display': 'inline-block'}),
                        html.Output(id='nextDose',children='0',style={'display': 'inline-block'}),
                    ]),
        html.Output(id='show_info',children=''),
        dcc.Graph(id='show_fig',figure={})
    ]),
],)
@app.callback(
    dependencies.Output('input-dose-container', 'children'),
    dependencies.Output('input-patients-container', 'children'),
    dependencies.Output('input-dlt-container', 'children'),
    dependencies.Input('add-filter', 'n_clicks'),
    dependencies.State('input-dose-container', 'children'),
    dependencies.State('input-patients-container', 'children'),
    dependencies.State('input-dlt-container', 'children'))
def display_input(n_clicks, children_dose,children_patients,children_dlt):
    new_input_dose = dcc.Input(
        id={
            'type': 'cohort_dose',
            'index': n_clicks
        }, value=1, type='number',min=1
    )
    children_dose.append(new_input_dose)
    new_input_patients = dcc.Input(
        id={
            'type': 'cohort_patient',
            'index': n_clicks
        },value=1,type='number'
    )
    children_patients.append(new_input_patients)
    new_input_dlt = dcc.Input(
        id={
            'type': 'cohort_dlt',
            'index': n_clicks
        },value=0,type='number'
    )
    children_dlt.append(new_input_dlt)
    return children_dose,children_patients,children_dlt


@app.callback(
    dependencies.Output('show_fig', 'figure'),
    dependencies.Input({'type': 'cohort_dose', 'index': dependencies.ALL}, 'value'),
    dependencies.Input({'type': 'cohort_patient', 'index': dependencies.ALL}, 'value'),
    dependencies.Input({'type': 'cohort_dlt', 'index': dependencies.ALL}, 'value'),
    dependencies.Input(component_id='n_doses', component_property='value'),
)
def display_output(values_dose,values_p,values_d,n_dose):
    # for i in
    dose = [value for (i,value) in enumerate(values_dose)]
    pat = [value for (i,value) in enumerate(values_p)]
    dlt= [value for (i,value) in enumerate(values_d)]
    n_cohort=len(pat)
    # print(pat)
    # print(dlt)
    # fig = go.Figure()
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    # n_dose = 5
    if len(pat)==0:
        cohort_max=1
    else:
        cohort_max= np.array(pat).max()
    # base_i =i*batch
    # base_i = [np.random.choice([i * 3 for i in range(n_dose)]) for i in range(n_cohort)]

    for i_dose in range(n_dose):
        fig.add_trace(go.Scatter(x=[0,n_cohort+3], y=[i_dose+1,i_dose+1],mode='lines',opacity=0,
                                 name='',),
                                 # name='Dose Level '+str(i_dose+1),),
                             # base=dose,
                             # marker_color='blue',
                      secondary_y=False,
                      )
    # fig.update_layout(showlegend=False,)
    # showlegend = False,
    fig.add_trace(go.Bar(x=[i+1 for i in range(n_cohort)], y=pat,
                         base=[i*cohort_max for i in dose],
                         marker_color='blue',
                         name='Positive response'),secondary_y=True,
                  )
    fig.add_trace(go.Bar(x=[i+1 for i in range(n_cohort)], y=dlt,
                         base=[i*cohort_max for i in dose],
                         marker_color='crimson',
                         name='Negative response'),secondary_y=True,
                  )
    fig.update_yaxes(range=[0, n_dose+1],secondary_y=False, title_text = "Dose Levels")
    # if len(pat)!=0:
    fig.update_yaxes(range=[0, (n_dose+1)*cohort_max],secondary_y=True,title_text = "Patients Number")
    fig.update_xaxes(tickvals=[i + 1 for i in range(len(pat))],range=[0,len(pat)+1],
                     title_text = "Cohort Number")

    fig.update_layout(barmode='stack', title_text='Dose-finding Observation',legend={'traceorder':'normal'})
    # text='test'+str(n_clicks)
    return fig

@app.callback(
    dependencies.Output(component_id='cMTD',component_property='children'),
    dependencies.Output(component_id='nextDose',component_property='children'),
    dependencies.Input(component_id='mc', component_property='n_clicks'),
    dependencies.State({'type': 'cohort_dose', 'index': dependencies.ALL}, 'value'),
    dependencies.State({'type': 'cohort_patient', 'index': dependencies.ALL}, 'value'),
    dependencies.State({'type': 'cohort_dlt', 'index': dependencies.ALL}, 'value'),
    dependencies.State(component_id='n_doses', component_property='value'),
    dependencies.State(component_id='target', component_property='value')
)
def update_fig(n_click,values_dose, values_p, values_d,num_doses,target):
    target_float = float(target)
    dose = [value for (i, value) in enumerate(values_dose)]
    pat = [value for (i, value) in enumerate(values_p)]
    dlt = [value for (i, value) in enumerate(values_d)]

    assign_info=np.zeros([num_doses])
    dlt_info=np.zeros([num_doses])
    for i in range(len(pat)):
        loc=dose[i]
        assign_info[loc]+=pat[i]
        dlt_info[loc]+=dlt[i]
    ndose = num_doses#len(values_p)

    tys = robjects.FloatVector(np.array(dlt))
    tns = robjects.FloatVector(np.array(pat))
    phi = target_float
    if len(values_p)>0:
        cidx=dose[-1]
        cMTD = robjects.r['abc_simu'](ndose,cidx,tys,tns,phi)
        cMTD = cMTD[0]
        if cMTD==99:
            nextDose='Early Stop'
        else:
            if cidx>cMTD:
                nextDose=cidx-1
            elif cidx==cMTD:
                nextDose=cidx
            else:
                nextDose=cidx+1
    else:
        cMTD='Please press the Run Simulation button'
        nextDose=' '

    return cMTD,nextDose

if __name__ == '__main__':
    app.run_server(debug=True,port=int(os.getenv('PORT', '4544')))