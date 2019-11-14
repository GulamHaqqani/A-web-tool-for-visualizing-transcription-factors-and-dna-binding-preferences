import io
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd
import numpy as np
from dash.dependencies import Input, Output, State
from scipy import stats
from itertools import product
import json
import base64

rows = ['STA','A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
cols = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','END']

data2 = pd.DataFrame(np.ones((21,21)), index=rows, columns=cols)

def all_tetra_subsets(ss=["A", "C", "G", "T"]):
    return [''.join(p) for p in product(ss, repeat=6)]

df_cols = all_tetra_subsets(["A", "C", "G", "T"])

colorscale = [
    [0.1,'#BA0A19'],
    [0.5,'#630109'],
    [1,'#F64A58']
]

colorscale1 = [
    [0, '#DDE3E3'],
    [0.1,'#FE5D53'],
    [0.5,'#FD3D30'],
    [0.75,'#EA1C0F'],
    [1,'#B41006']
]

#data = pd.read_csv("/Users/mac/Desktop/emissionmouseNOTTBOT2.csv", header=None)
#data.set_index(0,inplace=True)

#data.index = data.index.fillna('NA')
#data.columns = df_cols





app = dash.Dash(meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ])


app.layout = html.Div([
                html.Div([
                    html.H2("Visualizing Transcriptor Factors and DNA Binding Preferences", style={'margin-top':'7px'}), 
                   
                ], style={'backgroundColor':'grey','padding':'20px', 'width':'100%', 'height':'30px','text-align':'center', 'position':'fixed'}),
    
                html.Div([
                    
                    html.Div([
                    
                    html.Div([
                            
                            html.Div([
                                
                                dcc.Upload(id = "upload", children=html.Div(['Drag and Drop or ',
                                    html.A('Select File')
                                    
                                ]),style={
                                'width': '80%',
                                'height': '50px',
                                'lineHeight': '50px',
                                'borderWidth': '1px',
                                'borderStyle': 'dashed',
                                'borderRadius': '3px',
                                'textAlign': 'center',
                                'margin': '10px'
                            })
                                
                            ], style = {'margin-left': '15px', 'margin-top':'20px'}),
                        
                            html.Div([
                                html.Button('Upload File', id ='upload-file', style={'padding':'10px 22px', 'margin-top':'5px', 'margin-left':'15px', 'border-radius': '4px','font-size':'16px'})
                                
                            ], style = {'margin-left': '15px', 'margin-top':'5px'}),
                        
                        
                            html.Div(html.Label("DNA 6-mer" ),style = {'margin-top':'10px',
                                        'margin-left':'20px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                        'font-style': 'Playfair Display'}),
                            html.Div(
                             
                             dcc.Input(id = 'input_string', value = '', type='text', placeholder ='Please enter the DNA 6-mer',
                                      style ={
                                        'margin-top':'0px',
                                        'margin-left':'15px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                       'width': '200px',
                                       'border-radius': '4px',
                                       'font-size':'15px',
                                      }
                                      ), style = {'margin-left':'10px'}), 
                        
                             html.Div(html.Label("Protein Sequence" ),style = {'margin-top':'5px',
                                        'margin-left':'20px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                        'font-style': 'Playfair Display'}),
                            html.Div(
                             dcc.Input(id = 'input_seq', value = '', type='text', placeholder ='Please enter the protein sequence',
                                      style ={
                                        'margin-top':'0px',
                                        'margin-left':'15px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                       'width': '200px',
                                       'border-radius': '4px',
                                       'font-size':'15px',
                                      }
                                      ), style = {'margin-left':'10px'}), 
                            html.Div([
                            html.Button('Submit', id='button', style={'padding':'10px 26px', 'margin-top':'8px', 'margin-left':'15px', 'border-radius': '4px','font-size':'16px'}),
                            ], style = {'margin-left':'10px'}),
                            
                        
                             html.Div(html.Label("Select positive/negative 6-mers" ),style = {'margin-top':'5px',
                                        'margin-left':'20px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                        'font-style': 'Playfair Display'}),
                            html.Div([
                             dcc.Dropdown(id = 'dropdown', 
                                     options = [
                                         {'label': 'Positive 6-mers', 'value':'Pos'},
                                         {'label': 'Negative 6-mers', 'value':'Neg'},
                                     ], 
                                    searchable = False, 
                                    style = {
                                        'width':'200px',
                                        'padding':'8px 8px 5px',
                                        'margin-top': '-5px',
                                        'margin-left':'10px',
                                        'height':'30px',
                                        

                                    }
                                    )
                             ]),
                               html.Div(html.Label("Extreme values"),style = {'margin-top':'15px',
                                        'margin-left':'20px',
                                        'textalign' : 'left',
                                       'padding': '10px 8px',
                                        'font-style': 'Playfair Display'}),
                               html.Div([
                             dcc.Dropdown(id = 'MinMaxdropdown', 
                                     options = [
                                         {'label': 'Max value of a sequence', 'value':'Maxx'},
                                         {'label': 'Min value of a sequence', 'value':'Minn'},
                                     ], 
                                    searchable = False, 
                                    style = {
                                        'width':'200px',
                                        'padding':'8px 8px 5px',
                                        'margin-top': '-5px',
                                        'margin-left':'10px',
                                        'height':'30px',
                                        

                                    }
                                    )
                             ]),
                        
                            html.Div([
                                html.Button('Clear', id='clear', style={'padding':'10px 35px', 'margin-top':'8px', 'margin-left':'15px', 'border-radius': '4px','font-size':'16px'})
                            ], style = {'margin-left':'10px', 'margin-top':'20px'}),
                        
                        
                        
                            html.Div(id="Div", style = {'display':'none'}),
                            html.Div(id = 'Output2', style = {'margin-top':'50px','margin-left':'10px'})
                        
                            
                            
                       
            
                              
                             ], style = {'height':'630px',
                               'display':'inline-block',
                               'width':'25%',
                              'margin-left':'10px',
                              'margin-top':'5px', 
                              'position':'fixed',
                                'backgroundColor':'#F2F4F4'}),
                        
                    html.Div([dcc.Graph(id = 'map', style = {'width':'70%', 'display':'flex', 'margin-top':'3px'}), 
                             
                             
                             ],style = {'height':'630px',
                               'display':'inline-block',
                               'width':'75%',
                              'margin-left':'27%',
                              'margin-top':'5px', 
                              'position':'fixed',
                                'backgroundColor':'#F2F4F4'} )
                    ])
                    
                ], style = {'margin-top':'80px','backgroundColor':'#ABB2B9', 'width':'100%','height':'100%', 'position':'fixed'}),
    
])


def read_csv(contents, filename, date):
    content_type, content_string = contents.split(',')
    
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')), header=None)
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    return df




def max_min(radio, data):
    data_max_min = pd.DataFrame(np.zeros((21,21)),index=rows,columns=cols)
    k = 0
    for i in range(21):
        for j in range(21):
            if(radio == "Maxx"):
                data_max_min.iat[i,j] = max(data.iloc[k])
            else:
                data_max_min.iat[i,j] = min(data.iloc[k])
            k = k + 1

    hovertext = list()
    fig = go.Figure()
    for yi, yy in enumerate(cols):
        hovertext.append(list())
        for xi, xx in enumerate(rows):
            s  =   xx + yy
            max_s = max(data.loc[s])
            min_s = min(data.loc[s])
            Ser = data.loc[s]  #td[td==tmax].index[0]
            if(radio=="Maxx"):
                hovertext[-1].append('{}: {}'.format(Ser[Ser==max_s].index[0],np.around(max_s,5)))
                colorscale1 = [[0,'#FCB19B'],[0.25,'#FC5E2B'],[0.5,'#FA231D'],[0.75,'#D40A03'],[1,'#A60B06']]
            else:
                hovertext[-1].append('{}: {}'.format(Ser[Ser==min_s].index[0],np.around(min_s,5)))     
                colorscale1 = [[0,'#7295F7'],[0.5,'#3C6CF5'],[1,'#042B99']]
    
    figure = {'data':go.Heatmap(z=data_max_min.transpose(), x = rows, y = cols,colorscale=colorscale1, hoverinfo='text',text=hovertext,reversescale=False),
              'layout':go.Layout(height=600,width=900)
             }
    fig = figure
    return go.Figure(fig)
                
def graph_pos_neg(dropdown,data):
    hovertext = list()
    fig = go.Figure()
    for yi, yy in enumerate(cols):
        hovertext.append(list())
        for xi, xx in enumerate(rows):
            s  =   xx + yy
            df_values = data.loc[s].sort_values()
            df_index = df_values.index
            if(dropdown=="Neg"):
                color = [[0,"#2936DB"],[1,"#2936DB"]]
                hovertext[-1].append('{}: {}<br />{}: {}<br />{}: {}<br />{}: {}<br />{}: {}'.format(df_index[0],np.around(df_values[0],5),df_index[1],np.around(df_values[1],5),df_index[2],np.around(df_values[2],5),df_index[3],np.around(df_values[3],5),df_index[4],np.around(df_values[4],5)))        
            else:
                color = [[0,"#D51101"],[1,"#D51101"]]
                hovertext[-1].append('{}: {}<br />{}: {}<br />{}: {}<br />{}: {}<br />{}: {}'.format(df_index[-1],np.around(df_values[-1],5),df_index[-2],np.around(df_values[-2],5),df_index[-3],np.around(df_values[-3],5),df_index[-4],np.around(df_values[-4],5),df_index[-5],np.around(df_values[-5],5)))
    dt = pd.DataFrame(np.ones((21,21)), index=rows,columns=cols)
    figure = {'data':go.Heatmap(z=dt, x = rows, y = cols, xgap=0.5, ygap=0.5,hoverinfo='text',text = hovertext,colorscale=color, showscale=False), 
                 'layout':go.Layout(height=600, width = 900)}
    fig = figure
    return go.Figure(fig)
        
def graph_seq_highlighted(inp_str,inp_seq, data):
    dct = {'STA': 0,'A':1,'R':2,'N':3,'D':4,'C':5,'E':6,'Q':7,'G':8,'H':9,'I':10,'L':11,'K':12,'M':13,'F':14,'P':15,'S':16,'T':17,'W':18,'Y':19,'V':20, 'END':21}

    reversed_dct = dict(map(reversed, dct.items()))    
    fig = go.Figure()
    lst = list()
    spl_lst = list()
    
    spl_lst = list(inp_seq)
    lst.append('STA')
    [lst.append(s) for s in spl_lst]
    lst.append("END")

    fig = go.Figure()
    shapess = list()

    for i in range(len(lst) -1):
        t1 = dct[lst[i]],
        t0 = dct[lst[i+1]]
        if(t0 == 21 and lst[i+1] == "END"):
            yt = 19.55
            ytt = 20.5
        else:
            rr = reversed_dct[int(t0) + 1]
        temp = go.layout.Shape(
        type = "rect",
        x0 = t1[0] - 0.5,
        x1 = dct[rows[int(t1[0])]] - 0.5,
        y0 = yt if (lst[i+1]=="END") else dct[lst[i+1]] - 1.5,
        y1 = ytt if (lst[i+1]=="END") else dct[rr] - 1.5,
        line=dict(color = 'white', width=0.1),
        fillcolor='white'
        )
        shapess.append(temp)    
    
    if(inp_str!=""):
        df_2 = pd.DataFrame(np.zeros((21,21)), index=rows,columns=cols)
        
        data_3 = pd.DataFrame(data[inp_str])
        k = 0
        for i in range(21):
            for j in range(21):
                df_2.iat[i,j] = float(data_3.iloc[k])
                k = k + 1
        print("FROM WITHIN inp_STR location")
        #layout = fig.update_layout(shapes= shapess),
        if(inp_seq!=""):
            s = ""
            data_2 = pd.DataFrame(np.zeros((21,21)),index=rows, columns=cols)
            for i in range(len(lst)-1):
                s1 = lst[i]
                s2 = lst[i+1]
                num = df_2.loc[s1][s2]
                data_2.loc[lst[i]][lst[i+1]] = num
            fig = go.Figure(data = go.Heatmap(z = data_2.transpose(), x =rows,y=cols,colorscale = colorscale),layout=go.Layout(width=900, height=600) )
            #fig.update_layout(shapes=shapess)
            #go.Figure(fig)
    else:
        data_2 = pd.DataFrame(np.zeros((21,21)),index=rows, columns=cols)
        for i in range(len(lst)-1):
            s = lst[i] + lst[i+1]
            data_2.loc[lst[i]][lst[i+1]] = max(data.loc[s])
        fig = go.Figure(data = go.Heatmap(z=data_2.transpose(),x=rows, y=cols, colorscale=colorscale),layout=go.Layout(width=900,height=600))
        print(lst)
    return go.Figure(fig)        

@app.callback(Output('Div', 'children'),
              [Input('upload', 'contents')],
              [State('upload', 'filename'),
               State('upload', 'last_modified')])

def update_output(list_of_contents, list_of_names, list_of_dates):
    df_return = pd.DataFrame()
    if list_of_contents is not None:
        df_return = read_csv(list_of_contents, list_of_names, list_of_dates)
    return df_return.to_json()



@app.callback(Output('map','figure'),
             [Input('dropdown','value'),
              Input('button','n_clicks'),
             Input('MinMaxdropdown', 'value'),
              Input('Div', 'children'),
             Input('clear','n_clicks'),
             Input('upload-file','n_clicks')], 
             [State('input_string', 'value'),
             State('input_seq','value'),
             ])


def update_output(dropdown,click, radio,inter_val, clear, up, value, inp_seq):    
    fig = go.Figure()
    value = value.upper()
    inp_seq = inp_seq.upper()
    ret_seq_drawn = go.Figure()

    data = pd.DataFrame()
    if(up is not None):
        data = pd.read_json(inter_val)
        data.set_index(0,inplace=True)
        data.index = data.index.fillna('NA')
        data.columns = df_cols

    if(radio=="Maxx" or radio=="Minn"):
        fig = max_min(radio, data)
    elif(dropdown=="Pos" or dropdown=="Neg"):
        fig = graph_pos_neg(dropdown, data)
    elif((value=="" and inp_seq!="") or (value!="" and inp_seq!="")):
        fig = graph_seq_highlighted(value,inp_seq, data)
    elif(value==""):
        if(inp_seq!=""):
            fig = graph_seq_highlighted(value,inp_seq, data)
        else:
            dtt = pd.DataFrame(np.ones((21,21)), index=rows,columns=cols)
            figure = {'data':go.Heatmap(z=dtt, x = rows, y = cols, xgap=0.5, ygap=0.5,hoverinfo='x+y',colorscale="reds", showscale=False), 
                     'layout':go.Layout(height=600, width = 900)}
            fig = figure
    else:
        a = np.zeros(len(rows))
        df_2 = pd.DataFrame(a, index=rows)
        for i in range(21):
            df_2.insert(i, rows[i], a)
            
        df_2.drop(0,axis =1 ,inplace=True)
        data_3 = pd.DataFrame(data[value])
        
        k = 0
        for i in range(21):
            for j in range(21):
                df_2.iat[i,j] = float(data_3.iloc[k])
                k = k + 1
                
        figure = {'data':go.Heatmap(z = df_2.transpose(), x = rows, y = cols, colorscale=colorscale,reversescale=False),
                  'layout':go.Layout(height=600,width=900)}
        fig = figure

    return go.Figure(data=fig)


@app.callback([Output('dropdown','value'),
              Output('MinMaxdropdown','value'),
              Output('input_string', 'value'),
              Output('input_seq','value'),
              Output('button','n_clicks')],
             [Input('clear','n_clicks')])

def clear_update(n_clicks):
    return ("","","","","")

if __name__ == "__main__":
    app.run_server(debug=False)
