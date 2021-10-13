# coding=utf8

import os
import json
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_daq as daq
from dash import callback_context
import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.express as px
import pandas as pd
from flask import Flask
from dash import Dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dotenv import load_dotenv
#from .exceptions import ImproperlyConfigured
from collections import Counter
from numpy import random
from plotly.figure_factory import create_dendrogram
from plotly.subplots import make_subplots
import pickle
import bz2
import _pickle as cPickle
from pathlib import Path
import sqlite3
import urllib3

#TODO: Fix stderr  error message when changing the genome, reporting the missing column for hueing markers in the main graph

class ImproperlyConfigured(Exception):
    """Raise this exception when an environment variable is not set.
    """
    pass

####################### CONFIG ############################

path = Path(__file__).parent

DOTENV_PATH = os.path.join(os.path.dirname(__file__), ".env")
load_dotenv(DOTENV_PATH)

if not "DEBUG" in os.environ or "DYNO" in os.environ:
    # dyno if the app is on Heroku
    debug = False
# google analytics with the tracking ID for this app
# external_js.append('https://codepen.io/jackdbd/pen/rYmdLN.js')
else:
    debug = True
    dotenv_path = os.path.join(os.path.dirname(__file__), ".env")
    load_dotenv(dotenv_path)

try:
    py.sign_in(os.environ["PLOTLY_USERNAME"], os.environ["PLOTLY_API_KEY"])
except KeyError:
    raise ImproperlyConfigured("Plotly credentials not set in .env")


version_number = os.getenv("VERSION", "-1")

if not version_number:
    version_number = "1.0.0"

app_name = "Companion Dashboard - Acinetobacter Comparative Genomics - v" + version_number
server = Flask(app_name)

try:
    server.secret_key = os.environ["SECRET_KEY"]
except KeyError:
    raise ImproperlyConfigured("SECRET KEY not set in .env:")

external_js = []

external_stylesheets = [
    # dash stylesheet
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
    "https://fonts.googleapis.com/css?family=Lobster|Raleway",
    "//maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css", ]

app = Dash(name=app_name, server=server, external_stylesheets=[dbc.themes.LUMEN])

PLOTLY_LOGO = "https://applbio.biologie.uni-frankfurt.de/acinetobacter/wp-content/uploads/2017/11/for_logo.png"


#####################################################################

########################### DATA ####################################

http = urllib3.PoolManager()
resp = http.request('GET','https://aci-dash.s3.computational.bio.uni-giessen.de/data/p_feature_tables.pickle.bz2')
req_decompr = bz2.decompress(resp.data)
df = cPickle.loads(req_decompr)

df.index.rename("id", inplace=True)

relevant_columns = [  # "# feature",
   "genomic_accession",
   "assembly",
   "start",
   "end",
   "strand",
   "non-redundant_refseq",
   "name",
   "symbol",
   "locus_tag",
   "feature_interval_length",
   # "product_length",
   "attributes"]

df = df[(df["# feature"] == "CDS") & (df["class"] == "with_protein")][relevant_columns]

df.columns = [  # "feat",
   "Genomic Acc",
   "assembly",
   "Start",
   "End",
   "Str",
   "RefSeq Acc",
   "Protein Annotation",
   "Sym",
   "Locus_tag",
   "Len",
   # "product_len",
   "Comments"]

df['id'] = df.index

rel_path = 'https://aci-dash.s3.computational.bio.uni-giessen.de/data/extended_assembly2strain.csv'
genomes_df = pd.read_csv( rel_path , header="infer",
                         sep='\t', index_col=1, dtype=str)

genomes_dict = genomes_df.to_dict(orient="index")

rel_path = 'https://aci-dash.s3.computational.bio.uni-giessen.de/data/hogs2virulence_factors_with_source.tsv'
hog2vir_df = pd.read_csv(rel_path, header=None,
                         sep='\t',
                         index_col=0,
                         names=["query", "eval", "hit_id", "hit_description", "source"])

http = urllib3.PoolManager()
resp = http.request('GET','https://aci-dash.s3.computational.bio.uni-giessen.de/data/p_full_annot.pickle.bz2')
req_decompr = bz2.decompress(resp.data)
full_hog_table = cPickle.loads(req_decompr)

categories = ["baumannii(55)", "calcoaceticus(4)", "other_acb(34)", "haemolyticus(50)", "baylyi(9)",
          "lwoffii(71)", "brisouii(7)",
          "qingfengensis(4)"][::-1]

####################################################################

################### LAYOUT #########################################


def create_header():
    navbar = dbc.Navbar(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px"),
                                lg={'size': 1, 'offset': 0, 'order': 1},
                                ),
                        dbc.Col(dbc.NavbarBrand(app_name, className="ml-2"),
                                align="center",
                                style={'text-align':'center'},
                                lg={'size': 10, 'offset': 0, 'order': 2},),
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href="https://acinetobacter.de",
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
            # dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],
        # color="dark",
        # dark=True,
        style={'border': '1px solid rgba(0,0,0,0.125)'}
    )

    header = html.Header(navbar)
    return header

############# BODY ################

Card_genome_selection = dbc.Card([
        html.P(html.Strong("Search/Select Genome:"),
               style={"margin-bottom": "0.5rem",
                      "font-weight": 900}),
        dcc.Dropdown(
            id='genome-dropdown',
            options=[{'label': '{} - {}'.format(" ".join(v["species_name"].split('~', 2)[0:2]), k),
                      'value': k}
                     for k, v in genomes_dict.items()
                     ],
            value='GCF_000737145.1',
            clearable=False,
            searchable=True,
            optionHeight=35,
            placeholder="Select/Search Genome"
        ),
        html.Div(id="assembly-acc",
                 style={"display": "none"}
        ),
        dbc.Table([
            html.Tbody([
            html.Tr([html.Td(['NCBI Taxonomy ID:']), html.Td(id='taxid')]),
            html.Tr([html.Td(['Corrected Species:']), html.Td(id='assign')]),
            html.Tr([html.Td(['Sample Year:']), html.Td(id='year')]),
            html.Tr([html.Td(['Isolated from:']), html.Td(id='isol-site')]),
            html.Tr(
                [html.Td(['Publication:']), html.Td(html.A("Link", id='publication', target='_blank'))]),
            ]),
            ],
            style={"margin-top": "0.5rem"},
        ),
        ],
        body=True,
    )

Card_map = dbc.Card([
        #html.Label('Sampled:'),
        dcc.Graph(id="world-map",
                  # config={"modeBarButtonsToRemove": ['pan2d']}
                  config={"displayModeBar": False}
                  ),
        html.P(id="location-text")
        ],
        body=True,
    )

Card_data_table = dbc.Card(
                    dbc.CardBody([ html.Div([
                                        html.P("To filter protein entries, specifiy a string or number in the filter cell of each column. Text columns are searched for " + \
                                                   "entries containing the string. You may also use the following operators: eq (exactly equals), >, >=, <, and <= . " + \
                                                   "Text columns are compared by dictionary order. Numeric-only columns by value. " + \
                                                   "Use the checkboxes in the leftmost column to select one or more entries for the cards below.",
                                               style={"width":"75%",
                                                      "float":"left",
                                                      "margin-right": "20px",
                                                      #"border": "1px solid red"
                                                      },
                                               className="text-secondary"
                                               ),
                                        dbc.Button("Specify a list of accessions...", id="open-xl",
                                             # style={"width":"25%",
                                             #        "z-index":10,
                                             #        },
                                             className="btn btn-secondary",

                                                   ),
                                        ],
                                        style={"margin-bottom": 20,
                                                  #"float": "left",
                                                  "position": "relative",
                                                  "overflow": "auto",
                                                  "display": "inline-block",
                                                  }
                                    ),
                                dbc.Modal(
                                    [
                                        dbc.ModalHeader("Please specify a set of proteins accessions that you want to focus on:"),
                                        dbc.ModalBody(children=[
                                                                dbc.FormGroup(
                                                                    [
                                                                        #dbc.Label("Text"),
                                                                        dbc.Textarea(
                                                                            placeholder="WP_000446781.1, WP_001984992.1, WP_001020931.1, ...",
                                                                            id="user_protein_acc",
                                                                        ),

                                                                        dbc.FormText(
                                                                            "Comma or new-line separated list of protein RefSeq accessions of selected genome." + \
                                                                            " Not matching accessions will be ignored."
                                                                        ),
                                                                        html.Div(id='datatable_query_structure', style={'whitespace': 'pre', 'display':'none'})
                                                                    ]
                                                                )
                                                                ]
                                                      ),
                                        dbc.ModalFooter(
                                            dbc.Button("Show matches", id="close-xl", className="ml-auto")
                                        ),
                                    ],
                                    id="modal-xl",
                                    size="xl",
                                ),
                                dash_table.DataTable(
                                    id='datatable-interactivity',
                                    # columns=[
                                    #     {"name": i, "id": i, "deletable": True, "selectable": True} for i in df.columns
                                    #     if i != "assembly" if i != "id"
                                    # ],
                                    columns=[{"name": i, "id": i, "deletable": True, "selectable": True} for i in ["Genomic Acc","Start","End","Str","RefSeq Acc",
                                             "Protein Annotation","Sym","Locus_tag","Len","Comments"]
                                    ],
                                    #style_table={'overflowY': 'scroll'},
                                    # data=df.to_dict('records'),
                                    editable=False,
                                    filter_action="native",
                                    sort_action="native",
                                    sort_mode="multi",
                                    # column_selectable="single",
                                    row_selectable="multi",
                                    # row_deletable=True,
                                    # selected_rows=list(range(len(df.index))),
                                    page_action="native",
                                    page_current=0,
                                    page_size=15,
                                    style_cell_conditional=[
                                        {
                                            'if': {'column_id': "Sym"},
                                            #'textAlign': 'left'
                                            "width":"5%"},
                                        {
                                            'if': {'column_id': "Str"},
                                            'textAlign': 'center',
                                            "width": "5%"},
                                        {
                                        'if': {'column_id': "Protein Annotation"},
                                              # 'textAlign': 'left'
                                              "width": "30%"},
                                        {
                                        'if': {'column_id': "Start"},
                                              # 'textAlign': 'left'
                                              "width": "6%"},
                                        {
                                        'if': {'column_id': "End"},
                                              # 'textAlign': 'left'
                                              "width": "6%"},
                                        {
                                            'if': {'column_id': "Len"},
                                            # 'textAlign': 'left'
                                            "width": "6%"},
                                    ],
                                    style_data_conditional=[
                                        {
                                            'if': {'row_index': 'odd'},
                                            'backgroundColor': 'rgb(248, 248, 248)'
                                        }
                                    ],
                                    style_data={
                                        'whiteSpace': 'normal',
                                        'height': 'auto',
                                        'lineHeight': '18px'
                                    },
                                    style_cell={
                                        'overflow': 'hidden',
                                        'textOverflow': 'ellipsis',
                                        'maxWidth': 0,
                                        'font-size': '12px',
                                        'font-family': '"Source Sans Pro", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol"',

                                    },
                                    style_header={
                                        'backgroundColor': 'rgb(230, 230, 230)',
                                        'fontWeight': 'bold',
                                        'font-family': '"Source Sans Pro", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol"',
                                    },

                                    tooltip_header={'Len': 'Gene length',
                                                    'Locus_tag': 'Locus tag',
                                                    'Protein Annotation': 'Protein Annotation as in NCBI RefSeq',
                                                    'RefSeq Acc': 'NCBI RefSeq Protein Accession',
                                                    'Start': 'CDS feature Start position',
                                                    'End': 'CDS feature End position',
                                                    'Str': 'Relative strand',
                                                    'Genomic Acc': 'Contig/scaffold or replicon accession',
                                                    'Comments': 'Additional comments',
                                                    'Sym': 'Gene symbol'},

                                    # Style headers with a dotted underline to indicate a tooltip
                                    style_header_conditional=[{
                                        'if': {'column_id': col},
                                        'textDecoration': 'underline',
                                        'textDecorationStyle': 'dotted',
                                              } for col in ['Len', 'Locus_tag', 'Protein Annotation',
                                                            'RefSeq Acc', 'Start', 'End', 'Str',
                                                            'Genomic Acc', 'Comments', 'Sym']],
                                    #tooltip_delay=0,
                                    #tooltip_duration=None,
                                ),
                                html.Div(id='datatable-interactivity-container')
                              ],
                    ),
                    className="mb-3",
                )

Card_genome_browser = dbc.Card(
                        dbc.CardBody(
                            [
                                html.P("Coming Soon!", className="card-text"),
                                #dbc.Button("HO", color="danger"),
                            ],

                        ),
    className="mt-3",
    color = "light",
    outline = True,
)

Card_scatter_plot_and_controls = dbc.Card([
    html.Fieldset([
        html.Div([
            html.Div([
                html.Label('Taxonomic Range (Y-Axis):'),
                dcc.Dropdown(
                    id='y-axis',
                    options=[{'label': k,
                              'value': k}
                             for k in [cat for cat in full_hog_table.columns
                                       if '(' in cat
                                       ]
                             ],
                    value='complete_acb(93)',
                    clearable=False,
                    searchable=True,
                    optionHeight=25,
                    placeholder="Select Taxonomic Group"
                ),
            ],
            style={'width': '48%', 'display': 'inline-block'}
            ),
            html.Div([
                html.Label('Taxonomic Range (X-Axis):'),
                dcc.Dropdown(
                    id='x-axis',
                    options=[{'label': k,
                              'value': k}
                             for k in [cat for cat in full_hog_table.columns
                                       if '(' in cat
                                       ]
                             ],
                    value='other(141)',
                    clearable=False,
                    searchable=True,
                    optionHeight=25,
                    placeholder="Select Taxonomic Group"
                )
                ],
                style={'width': '48%', 'float': 'right', 'display': 'inline-block'}
            )
            ],
            className="form-group",
        ),

        html.Div([
            html.Label('Color By:',
                       style={'margin-right':'4px'}
            ),
            dcc.RadioItems(
                id='hue-criterion-radio',
                options=[
                    {'label': 'Genus-level Core/Pan Genome', 'value': 'aci_core231_of_234'},
                    {'label': "Phylostratum Gained", 'value': 'gained_at'},
                ],
                value='aci_core231_of_234',
                labelStyle={'display': 'inline-block', 'padding-left': '0.85rem',
                            'vertical-align': 'top'},
                style={'float':'right', 'display': 'block'},
                inputStyle={"margin-right": "0.3rem"}
            ),
            ],
            className="custom control custom-radio",
            style={#"min-height": "1.3125rem",
                    #"padding-left": "1.5rem",
                    #"box-sizing": "border-box",
                    "text-align": "left",
                    "display": 'inline-block'
                    #"line-height": 1.5
            }

        ),
        dcc.Checklist(
            options=[
                {'label': ' Highlight known virulence factors (hits in databases)', 'value': 'VIR'},
                # {'label': 'Montréal', 'value': 'MTL'},
                # {'label': 'San Francisco', 'value': 'SF'}
            ],
            id="highlights-checkb",
            value=[],
            inputStyle={"margin-right": "0.3rem"},
            labelStyle={'vertical-align': 'top',
                       }
        )
    ]
    ),
    html.Div(
        children=[
            dbc.Spinner([
                dcc.Graph(
                    id="graph-0",
                    config={"modeBarButtonsToRemove": ['toggleSpikelines', 'autoScale2d', 'hoverClosestCartesian']}
                )
            ], id="loading-spinner", color="primary", type="border"),  # Spinner
        ],
    ),
    html.Div([
        html.Div(
            daq.BooleanSwitch(
                id="jitter-option",
                on=False,
                label = {'label': "add jitter x ± [0, .25] to resolve overlaps",
                        "style": {'font-size':'12px'} #this must be part of a label object
                        },
                color='#158cba',
                labelPosition="right",
            ),
            className='custom-control custom-switch',
            style={'z-index': 1,
                   'width': '42%',
                   'display': 'inline-block',
                  }
        ),
        html.Div([
            html.Div(
                id='selected_data_points',
                style={"display": "none"}
            ),
            #dbc.Button('KUCKUCK')
            ],
            style={'z-index': 1,
                   'width': '42%',
                   'display': 'inline-block',
                   'float': 'right'}
        )
        ],
        style={'display':'inline-block', 'position':'relative', 'margin-top': '10px'}
    ),
    html.Div(html.P("Click on a datapoint to select a protein. On selection, enriched protein annotation and orthologous-group specific information " + \
             "are displayed in the cards shown below. Also, the selected protein and its genomic context are highlighted in the tabular view above. " + \
             "Beware, without jittering, a point may reflect several proteins sharing the exact same coordinates. " + \
             "Single-click on a legend entry deselect the class of data points, double-click focuses on it. " + \
             "On hover, you will find tools for navigation. Additionally, you may select a group of interest. Simply choose" + \
             " the 'lasso' or 'box select' tool, select the data points to open the filter window."),
             style={"width": "100%",
                    "float": "left",
                    "margin-right": "20px",
                    "margin-top": "20px",
                    "font-size": "0.875rem",
                    # "border": "1px solid red"
                    },
             className="text-secondary"
    ),
    ],
    body=True,
    style={'font-size': '12px'}
)


content_first_row = dbc.CardDeck(
    [Card_genome_selection,
     Card_map
     ],
    style = {"margin-bottom": 20,
             "margin-top": 20},
)

content_second_row = dbc.Row([
    dbc.Col([
        dbc.Tabs([
            dbc.Tab([
                Card_data_table,
            ],
                label="Proteins: Tabular View",
                label_style={"font-weight":900}
            ),
            dbc.Tab([
                Card_genome_browser,
            ],
                label="Proteins: Genome View",
            ),
        ]
        ),
    ],
        lg={'size': 12, 'offset': 0, 'order': 1}
    ),
],
    no_gutters=False,
)


Card_pan_accessory_pie_chart = dbc.Card([
                    html.Div(
                    [
                        #html.H1(children="Set Composition"),
                        dcc.Graph(id="set_composi_graph",
                                  style={"height": "100%"},
                                  config={"frameMargins": 0},
                                  figure = {
                                     "layout": {
                                         "title": "Composition",
                                         "height": 280,  # px
                                         "text-align": "center"
                                     },
                                 },
                        ),

                    ]
                    ),
                    ],
                    body=True,
                    style={"margin-bottom": 20},
                )

Card_sunburst_chart = dbc.Card([
                    html.Div(
                        [
                            # html.H1(children="Set Composition"),
                            dcc.Graph(id="set_composi_graph2",
                                      style={"height": "100%"},
                                      config={"modeBarButtonsToRemove": ['zoom2d',
                                                                         'zoomIn2d',
                                                                         'zoomOut2d',
                                                                         'hoverClosestCartesian',
                                                                         'resetViews',
                                                                         "autoScale2d",
                                                                         "resetScale2d"
                                                                         ]},
                                      figure={
                                      "layout": {
                                                    "title": "VIR Composition",
                                                    "height": 100,  # px
                                          "xaxis": {
                                              "visible": False
                                          },
                                          "yaxis": {
                                              "visible": False
                                          },
                                          "annotations": [
                                              {
                                                  "text": "Option not selected.",
                                                  "xref": "paper",
                                                  "yref": "paper",
                                                  "showarrow": False,
                                                  "font": {
                                                      "size": 12
                                                  }
                                              }
                                          ]
                                                },
                                      },
                                      ),

                        ]
                    ),
                ],
                    body=True
                )

Card_prevalence_bar_chart = dbc.Card([
                    dbc.CardHeader(
                        html.H4(
                            id="click-data",
                            ),
                    ),
                    dbc.CardBody(
                        html.Div(
                            [
                                #dcc.Graph(id="dendro",
                                #          config={"displayModeBar": False},
                                #          figure=dendro_fig
                                #),
                                dcc.Graph(id="protein-prevalence",
                                          config={"displayModeBar": False},
                                          figure={
                                              "layout": {
                                                            "title": "Clade Prevalences",
                                                            "height": 500,  # px
                                                  "xaxis": {
                                                      "visible": False
                                                  },
                                                  "yaxis": {
                                                      "visible": False
                                                  },
                                                  "annotations": [
                                                      {
                                                          "text": "Nothing selected. Click data point in graph.",
                                                          "xref": "paper",
                                                          "yref": "paper",
                                                          "showarrow": False,
                                                          "font": {
                                                              "size": 12
                                                          }
                                                      }
                                                  ]
                                                        },
                                              },
                                ),
                            ],
                            #style={'display': 'inline-block'},
                        ),
                    )
                    ]
                )

Card_additional_annotations_table = dbc.Card([
                    dbc.CardHeader(
                    html.H4(children="Enriched Protein Annotation"),
                    ),
                    dbc.CardBody([
                        html.Div(
                            [
                           # #html.H3(children="Set Composition"),
                           # dcc.Graph(id="genome_info_graph"),
                            ]
                        ),


                        dbc.Table([
                            html.Tbody([
                            html.Tr([html.Td(['NCBI RefSeq:']), html.Td(id='refseq',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['Predicted SCL class [GN]:']), html.Td(id='scl_pred',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['This class % across HOG:']), html.Td(id='scl_pc_across_hog',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['KEGG KO(s):']), html.Td(id='keggKO',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['KEGG Description(s):']), html.Td(id='kegg_description',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['COG ID:']), html.Td(id='cogid_1',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['COG Letter:']), html.Td(id='cog_letter_1',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['COG Description:']), html.Td(id='cog_description_1',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['VIR Factor Hit ID:']), html.Td(id='virulence_hit_patric_id',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['VIR evalue:']), html.Td(id='virulence_hit_evalue',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['VIR Description']), html.Td(id='virulence_hit_description',
                                                                        style={'padding-left':"5px"})],
                                    style={"border-bottom": "1px solid #ddd"}
                                    ),
                            html.Tr([html.Td(['VIR Source DB:']), html.Td(id='virulence_source',
                                                                        style={'padding-left':"5px"})]
                                    ),
                            ]),
                            ],
                            #bordered=True,
                            style={"margin-top": "0.5rem", "size":'lg'},
                            #responsive=True,
                            #striped=True,
                        ),
	                    ]
                    ),
                    html.Div(id='selected-data')
                    ],
                )
content_third_row =         dbc.Row([
            dbc.Col(
                Card_scatter_plot_and_controls,
                lg={'size': 9, 'offset': 0, 'order': 1},
                style={"margin-bottom": 20},
            ),
            dbc.Col([
                Card_pan_accessory_pie_chart,
                Card_sunburst_chart,
                ],
                lg={'size': 3, 'offset': 0, 'order': 2},
                # width=2
            ),
        ]
        )

content_fourth_row =         dbc.Row([
            dbc.Col(
                Card_prevalence_bar_chart,
                lg={'size': 6, 'offset': 0, 'order': 1}
            ),
            dbc.Col(
                Card_additional_annotations_table,
                lg={'size': 6, 'offset': 0, 'order': 2},
                            # width=2
            ),
            ],
            no_gutters=False,
        )

def create_content():
    content = html.Div(
        children=[
            content_first_row,
            content_second_row,
            content_third_row,
            content_fourth_row
        ],
    )
    return content


def create_footer():
    # footer_style = {"background-color": theme["background-color"], "padding": "0.5rem"}
    footer = dbc.Row(
        [
            dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px")),
            html.Section("DFG FOR2251 - Ebersberger Lab",
                         className="page-footer",
                         id='footer',
                         style={"margin": "20px"}),
        ],
        align="center")
    return footer


def serve_layout():
    layout = html.Div(
        children=[create_header(), create_content(), create_footer(),
                  ],
        className="container",
        style={'background-color': '#f2f2f2'}
    )
    return layout

app.layout = serve_layout
#server = app.server

# for js in external_js:
#    app.scripts.append_script({"external_url": js})
# for css in external_stylesheets:
#    app.css.append_css({"external_url": css})


################### STATIC ############################

CMAP = {'QI clade': px.colors.sequential.Greens[1],
        'BR clade': px.colors.sequential.Greens[2],
        'LW clade': px.colors.sequential.Greens[3],
        'BA clade': px.colors.sequential.Greens[4],
        'HA clade': px.colors.sequential.Greens[5],
        'ACB clade': px.colors.sequential.Greens[6],
        'BNS clade': px.colors.sequential.Greens[7],
        'B clade': px.colors.sequential.Greens[8],
        'Core': px.colors.sequential.Blues[3],
        'Accessory': px.colors.sequential.Blues[6],
        'Strain specific': "black",
        'other nodes': 'rgb(173,216,230)', #lightblue (soft cyan)
        'VIR': 'rgb(240, 173, 78)',
        'sampled_country': 'rgb(21,140,186)'  # blue
        }

@app.callback(
    [Output('world-map', "figure"),
     Output('taxid', 'children'),
     Output('assign', 'children'),
     Output('year', 'children'),
     Output('isol-site', "children"),
     Output('publication', 'href'),
     Output('location-text', "children"),
     Output('assembly-acc', "children")],
    [Input('genome-dropdown', 'value')])
def update_map(assembly_acc):
    # print("UPDATE_MAP", assembly_acc)
    gf = genomes_df.loc[assembly_acc]

    if pd.isnull(gf.loc["scope"]):
        gf.loc["scope"] = 'world'

    if pd.isnull(gf.loc["location"]):
        gf.loc["location"] = "unspecified place"

    if pd.isnull(gf.loc["publication"]):
        gf.loc["publication"] = '#'

    fig = px.choropleth(gf.to_frame().T, locations="iso_alpha3",
                        locationmode="ISO-3",
                        # color="lifeExp",  # lifeExp is a column of gapminder
                        # hover_name="country",  # column to add to hover information,
                        hover_data=["country", "location", "year"],
                        color = ['sampled_country'],
                        scope=gf.loc["scope"],
                        projection='equirectangular',
                        )

    fig.layout.update(showlegend=False,
                      height=280,
                      #title='Sampled in {}, {}'.format(gf.loc["location"], gf.loc["country"]),
                      margin=dict(l=0, r=0, t=0, b=0),
                      #config = {"frameMargins": 0}
                      )

    fig.update_traces(showscale=False)

    # fig.layout.update(height=700)
    return fig, gf.loc["tax_id"], gf.loc["additions_and_corrections"], \
           gf.loc["year"], gf.loc["sample_type"], gf.loc["publication"], \
           'Sampled at: {}, {}'.format(gf.loc["location"], gf.loc["country"]), \
           assembly_acc

@app.callback(
    #[
    [Output('datatable-interactivity', "data"),
     #Output('datatable-interactivity', "derived_virtual_selected_row_ids")],
    Output('graph-0', 'selectedData'),
     Output('datatable-interactivity','derived_filter_query_structure')],
    [Input('assembly-acc', "children")],
    [State("user_protein_acc", "value")])
def update_table(assembly_acc, selprots):
    # print("UPDATE_TABLE", assembly_acc, selprots )
    # con = sqlite3.connect("db_aci.sqlite")
    # sql = "SELECT * from feat_tables WHERE assembly = '{}';".format(assembly_acc)
    # dff = pd.read_sql_query(sql, con)
    dff = df[df["assembly"]==assembly_acc]

    if len(dff.index) < 1:
        return dff.to_dict("records"), None, None
            #, []
    return dff.to_dict("records"), None, None
        #, dff.index

@app.callback(
    [Output('set_composi_graph', 'figure'),
     Output('set_composi_graph2', 'figure'),
     Output('graph-0', 'figure')],
    [Input('hue-criterion-radio', "value"),
     Input('x-axis', "value"),
     Input('y-axis', "value"),
     Input('jitter-option', "on"),
     Input('highlights-checkb', "value"),
     Input('assembly-acc', "children"),
     Input('datatable-interactivity', 'derived_virtual_row_ids'),
     Input('datatable-interactivity', 'derived_virtual_selected_row_ids'),
     #Input('datatable-interactivity', 'active_cell')
     ])
def update_genome_info(hue_criterion, x_axis_category, y_axis_category, jitter, highlight,
                       assembly_acc, row_ids, selected_row_ids
                       #active_cell
                       ):
    # print("update_genome_info", hue_criterion, x_axis_category, y_axis_category, jitter, highlight, assembly_acc)
    # only selected genome
    dff = df[df["assembly"]==assembly_acc]
    #ctx = callback_context
    #print(ctx.triggered[0]['prop_id'])
    #    if ctx.triggered[0]['prop_id'] == 'assembly-acc.children':
    genome_set_size = len(dff.index)

    # only filtered rows
    if row_ids:
        dff = dff.loc[dff.index.intersection(row_ids)] # since .loc does not allow for missing indexes
    # only selected rows
    if selected_row_ids:
        dff = dff.loc[dff.index.intersection(selected_row_ids)]
    # if not active_cell is None:
    #     print("Active Cell:", active_cell)
    #     active_row_index = active_cell["row_id"]
    #     if not active_row_index in dff.index:
    #         active_row_index = None
    # else:
    #     active_row_index = None

    proteins = dff["RefSeq Acc"].tolist()
    hog_table = full_hog_table.loc[full_hog_table.index.intersection(proteins)]

    strain_specifics = set(proteins).difference(hog_table.index)
    strain_specific_vals = []
    for c in hog_table.columns:
        if c in ['gained_at', 'aci_core231_of_234', 'acb_core91_of_93']:
            strain_specific_vals.append("Strain specific")
        else:
            strain_specific_vals.append(0)
    ss_vals_array = [strain_specific_vals] * len(list(strain_specifics))
    hog_table = pd.concat([hog_table, pd.DataFrame(ss_vals_array,
                                                   columns=hog_table.columns,
                                                   index=list(strain_specifics))])

    # Pie Chart Prep

    stats = hog_table.groupby(hue_criterion).size().reset_index(name='counts')

    cmap = CMAP

    category_orders = {
        "aci_core231_of_234" : ['Core', 'Accessory','Strain specific']
    }

    if len(dff.index) == genome_set_size:
        pie_chart_header="#Proteins: {} (all)".format(len(dff.index))
    else:
        pie_chart_header="#Proteins: {}".format(len(dff.index))

    pie_fig = px.pie(stats, values='counts', names=hue_criterion, title=pie_chart_header,
                 color=hue_criterion, color_discrete_map=cmap)

    pie_fig.layout.update(
                      height=280,
                      #title='Sampled in {}, {}'.format(gf.loc["location"], gf.loc["country"]),
                      margin=dict(l=0, r=0, t=0, b=0),
                      legend=dict(
                            orientation="v",
                            yanchor="bottom",
                            y=-0.2,
                            xanchor="right",
                            x=1
                      ),
                      # showlegend=False,
                      )

    pie_fig.update_layout(
        font_family= '"Source Sans Pro", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol"',
        #font_size="1.3125rem",
        font_color='#222222',
        #title_font_family="Times New Roman",
        #title_font_color="red",
        #legend_title_font_color="green"
    )

    ## add jitter
    if jitter:
        hog_table[x_axis_category] = hog_table[x_axis_category].dropna().astype(int).apply(lambda n: n+(random.random_sample()-0.25))

    scatter_fig = px.scatter(hog_table,
                      y=y_axis_category, #"complete_acb(93)",
                      x=x_axis_category, #"other(141)",
                      color=hue_criterion,
                      marginal_y="box",
                      # marginal_x="box", #trendline="ols",
                      color_discrete_map=cmap,
                      custom_data=[hog_table.index],#needs to be as list as its arrayed(?)
                      template="simple_white"
                      )

    scatter_fig.layout.update(
        clickmode='event+select',
        #showlegend=False,
        #height=300,
        #title='Sampled in {}, {}'.format(gf.loc["location"], gf.loc["country"]),
        margin=dict(l=0, r=0, t=40, b=0),
        legend_title_text="",
        height=500
        #config = {"frameMargins": 0}
        )

    scatter_fig.update_traces(notched=False, selector=dict(type='box'))

    ### highlight vir factors
    if len(highlight) != 0: #list of possibly multiple features to be highlighted
        for i in highlight:
            a_list = [i if hog in hog2vir_df.index else float("NaN") for hog in hog_table['hog_id1'].to_list()]
            #a_list = {v: (i if v in hog2vir_df.index else nan) for v in hog_table['hog_id1'].to_list()}
            hog_table[i] = a_list
            hog_table[i].fillna(hog_table[hue_criterion], inplace=True)

        scatter_fig2 = px.scatter(hog_table,
                          y=y_axis_category,  # "complete_acb(93)",
                          x=x_axis_category,  # "other(141)",
                          color=i,
                          marginal_y="box",
                          # marginal_x="box", #trendline="ols",
                          color_discrete_map=cmap,
                          custom_data=[hog_table.index],  # needs to be as list as its arrayed(?)
                          template="simple_white"
                          )

        scatter_fig2.update_layout(clickmode='event+select')

        scatter_fig2.layout.update(
            margin=dict(l=0, r=0, t=40, b=0),
            legend_title_text="",
            height=500,
            # config = {"frameMargins": 0}
            )

        scatter_fig2.update_traces(notched=False, selector=dict(type='box'))

        stats2 = hog_table.groupby([hue_criterion, i]).size().reset_index(name='counts')

        sunburst_fig3 = px.sunburst(stats2, path=[hue_criterion,i], values='counts', color=i,
                           color_discrete_map=cmap)

        sunburst_fig3.layout.update(
            title= "VIR Composition",
            height=300,
            margin=dict(l=0, r=0, t=40, b=0),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.2,
                xanchor="right",
                x=1
            ),
            showlegend=True
        )

        return pie_fig, sunburst_fig3, scatter_fig2

            #hog_table[i] = hog_table_vir_subset = hog_table[hog_table["hog_id1"].isin(hog2vir_df.index)].index
        # colors = hog_table[hue_criterion].map(cmap)
        # hog_table_vir_subset = hog_table[hog_table["hog_id1"].isin(hog2vir_df.index)].index
        # #virs_cmap = {i:'rgb(240, 173, 78)' for i in hog_table_vir_subset}
        # #new_colors = colors.map(virs_cmap).fillna(colors)
        #
        # fig2.for_each_trace(
        #     lambda trace: trace.update(
        #         marker_color="rgb(240,173,78)") if trace.customdata in hog_table_vir_subset else (),
        # )
        # #fig2.for_each_trace( # accessory #core #strain specific are the traces
        # #    lambda trace: print(trace.__dict__),)
        #     #lambda trace: trace.update(marker_color=new_colors}
        # #print(new_colors)

    default_fig =  {
              "layout": {
                            "title": "VIR Composition",
                            "height": 100,  # px
                  "xaxis": {
                      "visible": False
                  },
                  "yaxis": {
                      "visible": False
                  },
                  "annotations": [
                      {
                          "text": "Option not selected.",
                          "xref": "paper",
                          "yref": "paper",
                          "showarrow": False,
                          "font": {
                              "size": 12
                          }
                      }
                  ]
                        },
              }

    return pie_fig, default_fig, scatter_fig


def create_normalized_barchart(series, identifier):

    if pd.isnull(series.loc["hog_id1"]):
        return {}

    categories = ["baumannii(55)", "calcoaceticus(4)", "other_acb(34)", "haemolyticus(50)", "baylyi(9)",
          "lwoffii(71)", "brisouii(7)",
          "qingfengensis(4)"][::-1]
    colors = ["#ff0000", "#e8e8ff", "#c8c8ff", "#ffd47f",'#00eeff', "#ff0088", "#4400ff", "#22ff00"][::-1]

    cmap={ele:colors[i] for i,ele in enumerate(categories)}

    total_counts = [55, 4, 34, 50, 9, 71, 7, 4][::-1]
    normalized = [round((int(series.loc[ele])/total_counts[i])*100,1) for i,ele in enumerate(categories)]
    data_dict = {identifier: normalized}
    data_df = pd.DataFrame.from_dict(data_dict, orient='index', columns=categories)
    data_df = data_df.T

    fig = px.bar(data_df, y=data_df.index, x=[identifier], title="Ortho-Family: " + series.loc["hog_id1"],
                orientation='h', text=data_df[identifier].to_list(), height=340,
                # color=data_df.index, color_discrete_map=cmap
                )

    #fig.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    fig.update_layout(uniformtext_minsize=8,
                      uniformtext_mode='hide',
                      margin=dict(l=0, r=0, t=40, b=0)
                      )

    fig.update_xaxes(title_text='Prevalence [%]',
                    # range=[0, 109]
                    )
    fig.update_yaxes(title_text='Taxonomic Group')

    fig.update_layout({
    'plot_bgcolor': 'rgba(0, 0, 0, 0)',
    'paper_bgcolor': 'rgba(0, 0, 0, 0)',
    })
    #fig.update_traces(marker_width=20)

    #fig.update_layout(showlegend=False)
    fig.layout.update(showlegend=False)

    return fig


@app.callback(
    [Output('click-data', 'children'),
     Output('protein-prevalence', 'figure'),
     Output('datatable-interactivity', 'active_cell'),
     Output('datatable-interactivity', "page_current"),
     Output('refseq', 'children'),
     Output('scl_pred', 'children'),
     Output('scl_pc_across_hog', 'children'),
     Output('keggKO', 'children'),
     Output('kegg_description', 'children'),
     Output('cogid_1', 'children'),
     Output('cog_letter_1', 'children'),
     Output('cog_description_1', 'children'),
     Output('virulence_hit_patric_id', 'children'),
     Output('virulence_hit_evalue', 'children'),
     Output('virulence_hit_description','children'),
     Output('virulence_source', 'children'),
    ],
    [Input('graph-0', 'clickData'),
     Input('assembly-acc', "children"),
     Input('datatable-interactivity', 'derived_virtual_row_ids')
     ])
def display_click_data(clicked_data, assembly_acc, row_ids):
    # print("DISPLAY_CLICKED_DATA", clicked_data, assembly_acc)
    # ctx = callback_context
    # print(ctx.triggered[0]['prop_id'])

    if clicked_data is None: #needs to be triggerd du to genome-dropdown (can be solved by storing as a DIV value)
        raise PreventUpdate()
        #return ["Compare Prevalences: No protein selected."], {}, None, 0, "No protein selected.", "", "", "", "", "", "", "", "", "", "", ""
    else:
        protein_acc = clicked_data["points"][0]["customdata"][0]
        #x = df.loc[protein_acc].loc["RefSeq Acc"]

        dff = df[df["assembly"]==assembly_acc]

        dff = dff.loc[dff.index.intersection(row_ids)]  # since .loc does not allow for missing indexes

        list_of_possible_row_indexes = dff[dff["RefSeq Acc"]== protein_acc].index
        if list(list_of_possible_row_indexes):
            active_row_index = list_of_possible_row_indexes[0] #select first if multiple
        else:
            raise PreventUpdate()
            #return ["Clade Prevalences: No protein selected."], {}, None, 0, "No protein selected.","","","","","","","","","","",""
        position = dff.index.get_loc(active_row_index)
        page = int(position/15)
        row = position % 15
        if row == 0 and page >=1:
            page -= 1

        try:
            hog_series = full_hog_table.loc[protein_acc] #f its not found it must be strain specific
        except:
            return ["Clade Prevalences: trivial for strain specific protein: " + protein_acc], \
                   {}, \
                   {'row': row, 'column': 4, 'column_id': 'RefSeq Acc', 'row_id': active_row_index}, \
                   page, \
                   "Strain specific protein selected.", "", "", "", "", "", "", "", "", "", "", ""

        this_series = full_hog_table.loc[protein_acc]
        this_series.fillna("N/A", inplace=True)

        ###postprocessing:

        if this_series['keggKO'] != 'N/A':
            this_series['keggKO'] = html.A(this_series['keggKO'],
                                           href="https://www.genome.jp/dbget-bin/www_bget?ko:{}".format(
                                               this_series['keggKO']),
                                           target='_blank')
        if this_series['cogid_1'] != 'N/A':
            this_series['cogid_1'] = html.A(this_series['cogid_1'],
                                            href="https://www.ncbi.nlm.nih.gov/research/cog/cog/{}".format(
                                                this_series['cogid_1']),
                                            target='_blank')

        this_series['refseq'] = html.A(protein_acc,
                                  href="https://www.ncbi.nlm.nih.gov/protein/{}".format(
                                      protein_acc),
                                target='_blank'
                                       )

        if this_series['virulence_hit_patric_id'] != 'N/A':
            this_series['virulence_hit_patric_id'] = html.A(this_series['virulence_hit_patric_id'][4:-1],
                                                         href='https://www.patricbrc.org/search/?keyword({})'.format(
                                                             this_series['virulence_hit_patric_id'][4:-1].replace('|',
                                                                                                                  '%7C')),
                                                        target='_blank'
                                                            )

        if this_series['virulence_hit_evalue'] != 'N/A':
            this_series['virulence_hit_evalue'] = format(float(this_series['virulence_hit_evalue']), '.2g')

        if this_series['scl_pc_across_hog'] != 'N/A':
            this_series['scl_pc_across_hog'] = round(this_series['scl_pc_across_hog'],2)

        #figure.addTraces(graphDiv, {
        #     x: X,
        #     y: Y,
        #     type: 'scatter',
        #     mode: 'markers',
        #     marker: {'color': 'black'},
        #     name: 'marker_trace'
        #   });
        #or
        # trace = go.Scattergl(
        #     x=x,
        #     y=y,
        #     mode='markers',
        #     marker=dict(
        #         size=10, )
        #     ,
        #     name='{}_{}'.format(plot_position, wtg),
        #     selected_marker_color='red')
        return ["Clade Prevalences for: " + protein_acc ], \
               create_normalized_barchart(hog_series, protein_acc), \
               {'row': row, 'column': 4 , 'column_id': 'RefSeq Acc', 'row_id':active_row_index}, \
                page, \
                this_series.loc['refseq'], \
                this_series.loc['scl_pred'].split("=")[-1], \
                this_series.loc['scl_pc_across_hog'], \
                this_series.loc['keggKO'], \
                this_series.loc['kegg_description'], \
                this_series.loc['cogid_1'], \
                this_series.loc['cog_letter_1'], \
                this_series.loc['cog_description_1'], \
                this_series.loc['virulence_hit_patric_id'], \
                this_series.loc['virulence_hit_evalue'], \
                this_series.loc['virulence_hit_description'], \
                this_series.loc['virulence_source']


def toggle_modal(n1, n2, is_open):
    # print("TOGGLE_MODAL",n1,n2,is_open)
    if n1 == -1:
        return False
    if n1 or n2:
        return not is_open
    return is_open

app.callback(
    Output("modal-xl", "is_open"),
    [Input("open-xl", "n_clicks"),
     Input("close-xl", "n_clicks")],
    [State("modal-xl", "is_open")],
)(toggle_modal)

@app.callback(
    [Output("user_protein_acc", "value"),
     Output("open-xl", "n_clicks"),
     Output("close-xl", "n_clicks")],
    [Input('selected_data_points', 'children')],
    State("close-xl", "n_clicks")
)
def trigger_modal(selected_data_points, close_n_clicks):
    # print("TRIGGER_MODAL", selected_data_points)
    if selected_data_points and len(selected_data_points) > 1:
        return selected_data_points, 1, close_n_clicks
    if selected_data_points == None or len(selected_data_points) == 1:
        raise PreventUpdate()
    return "", -1, 2

@app.callback(
### The problem here is, if we select data then the data gets filtered and the unselected are filtered from the graph
    Output('selected_data_points', 'children'),
    [Input('graph-0', 'selectedData')],
    State('selected_data_points', 'children')
)
def select_and_filter_selected_data(selectedData, selected_data_points):
    # print("SELECT_AND_FILTER", selectedData, selected_data_points)
    # ctx = callback_context
    # print(ctx.triggered[0]['prop_id'])
    # if ctx.triggered[0]['prop_id'] == 'assembly-acc.children':
    if selectedData == None or len(selectedData['points']) == 1 :
        if selected_data_points and not selected_data_points:
            return ""
        else:
            raise PreventUpdate()
    else:
        points = selectedData["points"]
        protein_accessions = [s["customdata"][0] for s in points]
        #dff = df[df["assembly"] == assembly_acc]
        #selected_row_ids = dff[dff["RefSeq Acc"].isin(protein_accessions)]["id"].to_list()
        #if len(selected_row_ids) == 1:
        #    PreventUpdate()
        #query = convert_list_to_filter_query(protein_accessions)
        return ", ".join(protein_accessions)
        #return selected_row_ids
        #return json.dumps(selectedData, indent=2)

# @app.callback(
#     [Output('graph-0', 'clickData'),
#     ],
#     [Input('datatable-interactivity', 'active_cell'),
#      Input('datatable-interactivity', 'derived_virtual_row_ids'),
#      Input('genome-dropdown', "value"),
#      Input('x-axis', "value"),
#      Input('y-axis', "value"),
#      ])
# def display_click_data(active_cell, row_ids, assembly_acc, xaxis, yaxis):
# #TODO: REVERSE select active cell that returns clickData as selected in the scatter plot, problem x, y coordinates
# #TODO: {'points': [{'curveNumber': 0, 'pointNumber': 2113, 'pointIndex': 2113, 'x': 122, 'y': 75, 'customdata': ['WP_001121112.1']}]}
# #TODO: problem is CIRCULAR DEPENDENCIES AND when jitter is added dynamically
#     dff = df[df["assembly"] == assembly_acc]
#     if row_ids:
#         dff = dff.loc[dff.index.intersection(row_ids)]  # since .loc does not allow for missing indexes
#
#     possible_pointNumbers = dff.get_loc(dff[dff["id"] == active_cell["points"][0]["customdata"][0]].index)
#     if list(possible_pointNumbers):
#         pointNumber = possible_pointNumbers[0]  # select first if multiple
#     else:
#         return None
#
#     return {'points': [{'curveNumber':0,
#                         'pointNumber':pointNumber,
#                         'pointIndex':pointNumber,
#                         'x': dff.iloc[pointNumber][xaxis],
#                         'y': dff.iloc[pointNumber][yaxis],
#                         'customdata': dff.iloc[pointNumber]["RefSeq Acc"]}]}

# @app.callback(
#     Output('selected-data', 'children'),
#     [Input('graph-0', 'selectedData')])
# def display_selected_data(selectedData):
#     return selectedData

# @app.callback(
#     Output('selected-data', 'children'),
#     [Input('basic-interactions', 'selectedData')])
# def display_selected_data(selectedData):
#     return json.dumps(selectedData, indent=2)


#IDEAD: multiindex the dataframe for assembly_accession if it takes too long

@app.callback(
    Output('datatable_query_structure', 'children'),
    [Input('close-xl', 'n_clicks')],
    [State('user_protein_acc', 'value')])
def update_output(n_clicks, value):
    # print("UPDATE_Output", n_clicks, value)
    if value == "" or not n_clicks:
        return None
    else:
        list_of_accessions = parse_input_accessions(value)
        query = convert_list_to_filter_query(list_of_accessions)
    return query

def parse_input_accessions(values_string):
    list_of_accessions = []
    lines = values_string.split('\n')
    for line in lines:
        [ list_of_accessions.append(e.strip()) for e in line.strip().split(',') ]
    return [l for l in list_of_accessions if l.strip()]

def convert_list_to_filter_query(values_list, column="RefSeq Acc"):
    # print("CONVERT", values_list)
    res = []
    for val in values_list:
        res.append('{{{0}}} eq {1}'.format(column, val))
    return ' or '.join(res) #gives an empty string if empty value_list

@app.callback(
    Output('datatable-interactivity', 'filter_query'),
    [Input('datatable_query_structure', 'children')]
)
def write_query(query):
    # print("WRITE_QUERY", query)
    if not query:
        return ''
    else:
        return query


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 9001))
    app.run_server(debug=debug, port=port, threaded=True)

