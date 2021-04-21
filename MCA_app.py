import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import dash.dependencies as dependencies
import os
import json
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# from mca_lib import get_simulation_samples,counts_
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# with open('./custom.geo.json', 'r') as f:
#     data = json.load(f)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

env_num=500
import numpy as np
def gen_env(num_doses,target_rate,env_num=1000):
    target_pos = range(num_doses)
    probs_2 = np.zeros([env_num*num_doses,num_doses],dtype='float32')
    for idx_ in range(env_num*num_doses):
        i= idx_%num_doses
        rand=[]
        for idx in range(num_doses):
            if idx<target_pos[i]:
                a=(target_rate-0.1)*np.random.rand()
            elif idx>target_pos[i]:
                a=(target_rate+0.1)+(1-target_rate-0.1)*np.random.rand()
            else:
                a=target_rate
            rand.append(a)
        rand_ = np.array(rand)
        probs_2[idx_,:] = rand_
    probs = probs_2
    probs = np.sort(probs,axis=-1)
    return probs
def gen_env_fig(num_doses,env,indices=None,flag=True):
    # marker_size=5
    # num_doses = 5
    titles=()
    # titles = ('Event Type '+ str(i+1) for i in range(num_doses))
    for i in range(num_doses):
        titles+=('Event '+ str(i+1),)
    if not flag:
        fig = make_subplots(rows=num_doses, cols=1,
                            shared_xaxes=True,
                            vertical_spacing=0.05,subplot_titles=titles)
        colors_list = plotly.colors.DEFAULT_PLOTLY_COLORS[:num_doses]
        opacity = 1  # 0.5
        marker_size = 5
        sample_num = env_num
        for row_i in range(num_doses):
            x_data_i = env[row_i::num_doses]
            # if row_i < 4:
            #     sort_index = np.argsort(x_data_i[:, row_i + 1])
            # else:
            #     sort_index = np.argsort(x_data_i[:, row_i - 1])
            # print(sort_index[:4])
            if row_i != 1:
                sort_index = np.argsort(x_data_i[:, 1])
            else:
                sort_index = np.argsort(x_data_i[:, 2])
            # if row_i == 0 or row_i == 1:
            # sort_index = np.argsort(x_data_i[:, 1])
            x_data_i = x_data_i[sort_index]
            for i in range(sample_num):
                data_i = x_data_i[i]
                for i_dose in range(num_doses - 1):
                    fig.add_trace(
                        go.Scatter(x=[data_i[i_dose], data_i[i_dose + 1], ], y=[i for j in range(2)], mode='lines',
                                   marker_color=colors_list[i_dose], line_width=1,
                                   opacity=opacity), row=row_i + 1, col=1, )
            fig.update_yaxes(title_text='Scenarios', row=row_i + 1, col=1)
        fig.update_xaxes(title_text="Toxicity Probabilities for Different Dose levels",row=num_doses, col=1,)
        fig.update_layout(showlegend=False,title_text="Generated Random Scenarios for Monte Carlo Simulation")
        fig.update_layout(height=800)
        # fig.update_layout(width=int(width),height=)
        return fig
    else:
        # print(indices.size())
        # indices_=indices.view(-1)
        indices_ = indices.reshape(-1)
        # indices_=list(indices_.numpy())
        indices_=list(indices_)
        print(len(indices_))
        flag=0.05*np.ones([env_num*num_doses])
        flag[indices_]=1
        fig = make_subplots(rows=num_doses, cols=1,
                            shared_xaxes=True,
                            vertical_spacing=0.05,subplot_titles=titles)
        colors_list = plotly.colors.DEFAULT_PLOTLY_COLORS[:num_doses]
        opacity = 1  # 0.5
        marker_size = 5
        sample_num = env_num
        for row_i in range(num_doses):
            x_data_i = env[row_i::num_doses]
            flag_i = flag[row_i::num_doses]
            if row_i != 1:
                sort_index = np.argsort(x_data_i[:, 1])
            else:
                sort_index = np.argsort(x_data_i[:, 2])
            x_data_i = x_data_i[sort_index]
            flag_i = flag_i[sort_index]
            for i in range(sample_num):
                data_i = x_data_i[i]
                flag_ii=flag_i[i]
                for i_dose in range(num_doses - 1):
                    fig.add_trace(
                        go.Scatter(x=[data_i[i_dose], data_i[i_dose + 1], ], y=[i for j in range(2)], mode='lines',
                                   marker_color=colors_list[i_dose], line_width=1,
                                   opacity=flag_ii,), row=row_i + 1, col=1, )
            fig.update_yaxes(title_text='Scenarios', row=row_i + 1, col=1)
        fig.update_xaxes(title_text="Toxicity Probabilities for Different Dose levels", row=num_doses, col=1, )
        fig.update_layout(showlegend=False, title_text="Generated Random Scenarios for Monte Carlo Simulation")
        fig.update_layout(height=800)
        return fig
def get_simulation_samples(num_doses,trial_num, simu_envs):
    # trial_num=torch.from_numpy(trial_num)
    # simu_envs=torch.from_numpy(simu_envs)
    for i in range(num_doses):
        probs = simu_envs[:, i]
        # print(trial_num)
        # print(probs)
        trial_num_tmp=trial_num[i]*np.ones_like(probs)
        trial_num_tmp = trial_num_tmp. astype(int)
        samples = np.random.binomial(n=trial_num_tmp,p=probs)
        # bern = Binomial(total_count=trial_num[i], probs=probs)
        # samples = bern.sample()
        samples = samples.reshape(-1,1)
        # samples = samples.view(-1, 1)
        if i == 0:
            samples_all = samples
        else:
            samples_all = np.concatenate([samples_all, samples],axis=-1)
            # samples_all = torch.cat([samples_all, samples], dim=-1)
        # samples_all=torch.from_numpy(samples_all)
    return samples_all
def counts_(num_doses,samples, possible_state):
    # possible_state = torch.from_numpy(possible_state)
    # counts = torch.zeros([num_doses, ])
    counts = np.zeros(([num_doses,]))
    diff = samples - possible_state[None, :]
    # np.abs()
    diff = np.abs(diff).sum(-1)
    # print(diff[:10],diff.size())
    indices = (diff == 0).nonzero()
    # print(len(indices))
    # print(indices)
    indices = indices[0]
    # print(indices)
    # dd
    # TODO:env sensitive
    for idx in indices:
        count_idx = idx % num_doses
        counts[count_idx] = counts[count_idx] + 1
    return counts,indices

app.layout = html.Div(children=[
    html.H1(children='Monte Carlo Simulation for Dose-finding'),
    html.Div(
        children=[
            html.Label('Target:',style={'display': 'inline-block'}),
            dcc.Input(id='target', value='0.3',type='text',style={'display': 'inline-block'}),
            # dcc.Input(0.3,id='target'),
            html.Label('Number of dose levels:',style={'display': 'inline-block'}),
            dcc.Input(id='n_doses',value=5,type='number',min=2,style={'display': 'inline-block'}),
            # html.Button(children='Generate Monte Carlo Environments',id='gen_button'),
            # dcc.Graph(id='envs',figure=fig_),
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
        html.Button(children='Run MC Simulation',id='mc',n_clicks=0),
        html.Output(id='show_info',children=''),
        dcc.Graph(id='show_fig',figure={})
        # html.Div(id='dropdown-container-output')
    ]),
    html.Div(children=[dcc.Graph(id='fig-bar',figure={}),
                        dcc.Graph(id='fig-selected',figure={})
    ])
],)

# actual_trial_num = torch.zeros([num_doses, ])
# state_test = torch.zeros([num_doses, ])
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
        # options=[{'label': i, 'value': i} for i in ['NYC', 'MTL', 'LA', 'TOKYO']]
    )
    children_dose.append(new_input_dose)
    new_input_patients = dcc.Input(
        id={
            'type': 'cohort_patient',
            'index': n_clicks
        },value=1,type='number'
        # options=[{'label': i, 'value': i} for i in ['NYC', 'MTL', 'LA', 'TOKYO']]
    )
    children_patients.append(new_input_patients)
    new_input_dlt = dcc.Input(
        id={
            'type': 'cohort_dlt',
            'index': n_clicks
        },value=0,type='number'
        # options=[{'label': i, 'value': i} for i in ['NYC', 'MTL', 'LA', 'TOKYO']]
    )
    children_dlt.append(new_input_dlt)
    return children_dose,children_patients,children_dlt


@app.callback(
    # dependencies.Output('show_info', 'children'),
    dependencies.Output('show_fig', 'figure'),
    dependencies.Input({'type': 'cohort_dose', 'index': dependencies.ALL}, 'value'),
    dependencies.Input({'type': 'cohort_patient', 'index': dependencies.ALL}, 'value'),
    dependencies.Input({'type': 'cohort_dlt', 'index': dependencies.ALL}, 'value'),
    dependencies.Input(component_id='n_doses', component_property='value'),
    # dependencies.Input('show', 'n_clicks'),
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
    # dependencies.Output(component_id='output',component_property='children'),
    dependencies.Output(component_id='fig-bar',component_property='figure'),
    dependencies.Output(component_id='fig-selected',component_property='figure'),
    dependencies.Input(component_id='mc', component_property='n_clicks'),
    dependencies.State({'type': 'cohort_dose', 'index': dependencies.ALL}, 'value'),
    dependencies.State({'type': 'cohort_patient', 'index': dependencies.ALL}, 'value'),
    dependencies.State({'type': 'cohort_dlt', 'index': dependencies.ALL}, 'value'),
    dependencies.State(component_id='n_doses', component_property='value'),
    dependencies.State(component_id='target', component_property='value')
)
def update_fig(n_click,values_dose, values_p, values_d,num_doses,target):
    # dff = df[df.year==year]
    # fig_geo = gen_fig(dff, data, column)
    target_float = float(target)
    env_ = gen_env(num_doses, target_float, env_num=500)
    # fig_ = gen_env_fig(num_doses, env_)
    dose = [value for (i, value) in enumerate(values_dose)]
    pat = [value for (i, value) in enumerate(values_p)]
    dlt = [value for (i, value) in enumerate(values_d)]

    assign_info=np.zeros([num_doses])
    dlt_info=np.zeros([num_doses])
    for i in range(len(pat)):
        loc=dose[i]
        assign_info[loc]+=pat[i]
        dlt_info[loc]+=dlt[i]
    # assign_info=get_obs_info(assign)
    # dlt_info=get_obs_info(dlt)
    # print(assign_info)
    # print(dlt_info)
    samples = get_simulation_samples(num_doses, assign_info, env_)
    counts, indices = counts_(num_doses, samples, dlt_info)
    if assign_info.sum()==0:
        fig_selected = gen_env_fig(num_doses, env_, flag=False)
    else:
        fig_selected = gen_env_fig(num_doses, env_, indices=indices)
    # print(counts)
    # print(indices[:10])
    # actual_trial_num = assign
    # output='test'#+assign_info[0]
    fig_bar=px.bar(x=np.arange(1,num_doses+1),y=counts,labels={'x':'Event Type','y':'Number of Occurrences'})
    fig_bar.update_xaxes(tickvals=[i + 1 for i in range(num_doses)],)

    return fig_bar,fig_selected


def get_obs_info(text):
    assign_info = text.split(',')
    ret = np.zeros([len(assign_info)])
    for i in range(len(assign_info)):
        item = assign_info[i].split(':')
        ret[i] = int(item[1])
    return ret

if __name__ == '__main__':
    app.run_server(debug=True,port=int(os.getenv('PORT', '4544')))