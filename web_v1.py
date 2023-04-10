import streamlit as st
import pandas as pd
import plotly.express as px
import requests


#Titles
st.title('Caldendrin binding partners analysis')
st.header('Swetha Iyer, Lisa Phan')
st.subheader('Bioinformatics Spring 2023')

######### Links ########

#make sure that raw = true
link = r'https://github.com/Lisa-Phan/Bioinfo_Project_Spring2023/blob/master/movie2.mp4?raw=true'
video_file = requests.get(link).content

# Display the video in the Streamlit app
st.video(video_file)

def introduction():
    st.subheader('The role of Caldendrin')


########### Plot1 ###########
def plot_interpro():
    st.title('InterPro domains of EF-hand 1 and 2 binders')
    data = pd.read_csv('Interpro_binder_annotation.tsv', sep = '\t')
    data = data.sort_values(by = 'n', ascending = False)
    n_slider = st.slider('Select the minimum value of n:', min_value=0, max_value=int(data['n'].max()), value=15)

    # Filter data based on slider value
    filtered_data = data.query("n >= @n_slider")

    fig = px.bar(
        filtered_data,
        x="Interpro.Description",
        y="n",
        title='Enriched InterPro Domains of EF-hand 1 and 2 binders',
        hover_data=['InterPro'],
        color="n"
    )

    fig.update_xaxes(tickangle=45, title_text='Interpro Description', title_font={"size": 20}, title_standoff=25)
    fig.update_layout(width=1200,height=800, margin=dict(t=0, b=0, l=0, r=0))

    st.plotly_chart(fig, theme='streamlit', use_container_width=True)


########### Plot2 ###########
#Biomart query

    st.title('InterPro domains of post-synaptic density proteins')
    Interpro_PSD_link = r"psd_interpro_count.tsv"
    psd = pd.read_csv(Interpro_PSD_link, sep = '\t')
    psd = psd.sort_values(by = 'n', ascending = False)

    #rename columns
    psd = psd.rename(columns = {'Interpro Description':'Interpro_Description', 
                                'Interpro ID':'Interpro_ID',
                                'n':'n_PSD', 
                                'index': 'InterPro'})

    n_slider_2 = st.slider('Select the minimum value of n:', min_value=0, max_value=int(psd['n_PSD'].max()), value=15)
    psd_filtered = psd.query("n_PSD >= @n_slider_2")


    fig2 = px.bar(
        psd_filtered,
        x="Interpro_Description",
        y="n_PSD",
        title='InterPro Domains of Post-synaptic density proteins',
        hover_data=['InterPro'],
        color="n_PSD"
    )
    fig2.update_xaxes(tickangle=45, title_text='Interpro Description', title_font={"size": 20}, title_standoff=25)
    fig2.update_layout(width=1200,height=800, margin=dict(t=0, b=0, l=0, r=0))
    st.plotly_chart(fig2, theme='streamlit', use_container_width=True)

########### Plot3 ###########
#correlation plot


    st.title('Proportion of Interpro domains \ncomparing EF-hand 1 & 2 binders vs. total post synaptic density proteins\n')
    #calculate the proportion of each InterPro domain in the binder set
    data['prop'] = data['n']/data['n'].sum()

    #calculate the proportion of each InterPro domain in the PSD set
    psd['prop'] = psd['n_PSD']/psd['n_PSD'].sum()

    #merge the two dataframes
    prob = pd.merge(data, psd, on = 'InterPro', how = 'outer')


    search_term = st.text_input('Search for a data point from IPR number:')
    if search_term:
        prob = prob.loc[prob['InterPro'] == search_term]

    fig3 = px.scatter(prob, x = 'prop_x', 
                    y = 'prop_y', 
                    hover_data = ['Interpro.Description', 'Interpro_ID'], 
                    color = 'prop_x')

    fig3.update_layout(width=1200,height=800, margin=dict(t=0, b=0, l=0, r=0))
    fig3.update_xaxes(title_text='Proportion in EF-hand 1 & 2 binders)', title_font={"size": 20}, title_standoff=25)
    fig3.update_yaxes(title_text='Proportion in PSD proteins', title_font={"size": 20}, title_standoff=25)
    st.plotly_chart(fig3, theme='streamlit', use_container_width=True)

    #Interactors
    B_interactors_link = r"result_2.csv"
    B_interactors = pd.read_csv(B_interactors_link)

    #Create a barplot from the ID columns count

    # subset the columns of interest
    B_int = B_interactors[['#Query', 'ID']]
    #print(new_df)
    ID_count = B_int['ID'].value_counts().reset_index().drop([0]).rename(columns = {'index':'Domains', 'ID':'Count'})
    

    #interactive slider
    n_slider_3 = st.slider('Select the minimum value of n:', min_value=0, max_value=int(ID_count['Count'].max()), value=15)
    ID_count_fil = ID_count.query("Count >= @n_slider_3")

    #interactive search
    search_term_2 = st.text_input('Search for domains based on name:')
    if search_term_2:
        ID_count_fil = ID_count[ID_count.Domains.str.contains(search_term_2, regex= True, na=False)]
    #foo[foo.b.str.contains('oo', regex= True, na=False)]
    fig4 = px.bar(
        ID_count_fil,
        x='Domains',
        y='Count',
        title='InterPro Domains of Post-synaptic density proteins',
        hover_data=['Domains'],
        color="Count"
    )
    st.plotly_chart(fig4, theme='streamlit', use_container_width=True)

# Define dropdown menu options
options = ['Select an option', 'InterPro domains distribution']

# Create dropdown menu
choice = st.sidebar.selectbox('Choose an option:', options)

# Display selected option
if choice == 'Introduction':
    introduction()
elif choice == 'InterPro domains distribution':
    # Create slider for minimum value of 'n'
    # Create plot
    plot_interpro()





