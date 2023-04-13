import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import plotly.express as px
import requests

########### Global functions ############
#########################################

def color_marker(known_binding, interpro_col):
    """ 
    Param: known_binding_interpro: a list of ID to color against, interactor_df['InterPro']
    Param: interpro_col: column to check against
    Return: a list of colors same length as interpro_col
    """   
    return ['Known_binder' if x in known_binding else 'Other' for x in interpro_col]

########## End of global functions #######
##########################################


##############################################
######### Global plotting variables ##########

#data is all proteins that binds EF-hand 1 and 2
binder_link = r'Interpro_binder_annotation.tsv'
data = pd.read_csv(binder_link, sep = '\t')
data = data.sort_values(by = 'n', ascending = False)

Interpro_PSD_link = r"psd_interpro_count.tsv"
#psd is proteinst that exists in the PSD
psd = pd.read_csv(Interpro_PSD_link, sep = '\t')
psd = psd.sort_values(by = 'n', ascending = False)

#rename columns
psd = psd.rename(columns = {'Interpro Description':'Interpro_Description', 
                            'Interpro ID':'Interpro_ID',
                            'n':'n_PSD', 
                            'index': 'InterPro'})

data['prop'] = data['n']/data['n'].sum()
psd['prop'] = psd['n_PSD']/psd['n_PSD'].sum()

#merge the two dataframes
prob = pd.merge(data, psd, on = 'InterPro', how = 'inner')

B_interactors_link = r"result_2.csv"
B_interactors = pd.read_csv(B_interactors_link)
ID_count = B_interactors['ID'].value_counts().reset_index().drop([0]).rename(columns = {'index':'Domains', 'ID':'Count'})

interactor_link = r'known_caldendrin_binding_protein.txt'
interactor_df = pd.read_csv(interactor_link, sep = '\t')

#Modifying the interactor data
prob['color_marker'] = color_marker(interpro_col = prob['InterPro'].to_list(), known_binding = interactor_df['InterPro'].to_list())

binder_dict = dict(zip(interactor_df['InterPro'].to_list(), interactor_df['Protein'].to_list()))
prob['Binder_Protein_name'] = prob['InterPro'].map(binder_dict, na_action = 'ignore')

######## End of global plotting variables ########
##################################################




################ Assets Links ################
#############################################
link = r'https://github.com/Lisa-Phan/Bioinfo_Project_Spring2023/blob/master/movie2.mp4?raw=true'
image_venn = r'venn.png'
image_scheme = r'schematic.png'
image_motif = r'Meme_motif_binder_protein.png'


################# End of links ###############
#############################################

def intro():
    #Titles
    st.title('Identifying Caldendrin binding partners')
    st.header('Swetha Iyer, Lisa Phan')
    st.subheader('Bioinformatics Spring 2023')

    video_file = requests.get(link).content

    # Display the video in the Streamlit app
    st.video(video_file)

    #Description text
    st.subheader('Introduction')
    intro_cont = st.container()
    intro_cont.write("**Caldendrin** is a member of the neuronal calcium-binding protein (nCABP1) \
                     family that is crucial for the growth and development of neurons in the \
                     peripheral nervous system. It is known to interact with voltage-gated calcium channels (CaVs) \
                     and modulate their activity (Tippens et al., 2007). Caldendrin is particularly enriched \
                     in the postsynaptic density (PSD), a large protein complex involving more than 1000 proteins \
                     (Kreutz et al., 2019; Takumi et al., 2018). Despite being a highly ubiquitous protein, \
                     the interactions between Caldendrin and other PSD components are poorly understood.")

    # intro_cont.write("In this project, we seek to perform in silico protein interaction analysis to identify \
    #                  potential partners of Caldendrin in the vast cloud of proteins found in the PSD. \
    #                  Assuming that Caldendrin and other proteins within the calcium binding protein \
    #                  family possess common motif preferences, discovery of these motifs would allow us \
    #                  to better understand Caldendrin’s biological network .")

    # intro_cont.write('Using FoldSeek and the AlphaFold database,\
    #                  we will first generate a series of homologous models that can \
    #                  be used as a lead whose known binding partners are potential ligands.\
    #                  Using annotations from BioGRID and literature search, along with subcellular \
    #                  localization and expression profile filters, we will compile a list of these \
    #                  ligand proteins. Next, we will identify the motifs in this dataset using \
    #                  readily available online tools. To assess whether the proposed pipeline \
    #                  can identify Caldendrin binding partners, we will compare the motif profiles \
    #                  with those found in sequences of 7 known Caldendrin binding proteins. \
    #                  As an additional benchmark, we will compare our results with the output of ISLAND, \
    #                  a machine learning sequence-based protein binding affinity predictor (Andleeb & Minhas et al., 2020)')

    intro_cont.write("In this project, we seek to perform in silico protein interaction analysis to identify \
                     potential partners of Caldendrin in the vast cloud of proteins found in the PS. \
                     Since motifs-motifs interactions are known to \
                     re-occur across many distinct protein pairs, we hypothesize that matching Caldendrin's binder motif (**binder A**) \
                     to their ligand motifs (**ligand B**) will reveal clues of what proteins are likely to be Caldendrin's partners. ")

    st.image(image_scheme, use_column_width=True)
    intro_cont.write("In Caldendrin, the set of motifs in **binder A** consists of a long disordered N-terminus (not present in PDB above) \
                     and a series of EF-hand domains. \
                     In order to obtain the **ligand B** set, we will query for proteins that are known binders of **motif A** through InterPro and \
                     Pfam annotations offered by Uniprot. We will limit our total pool of potential candidates \
                     based on GO annotations for subcellular localization in PSD. To assess whether the proposed pipeline \
                     can identify Caldendrin binding partners, we will compare the motif profiles with \
                     those found in sequences of 7 known Caldendrin binding proteins. ")
    
    st.image(image_venn, use_column_width=True)

    citation = st.container()
    citation.write("**1.** Mundhenk, J., Fusi, C., & Kreutz, M. R. (2019). \
                   Caldendrin and calneurons—EF-hand CaM-like calcium \
                   sensors with unique features and specialized neuronal \
                   functions. Frontiers in Molecular Neuroscience, 12, 16.")
    citation.write("**2.** Kaizuka, T., & Takumi, T. (2018). Postsynaptic density  \
                   proteins and their involvement in neurodevelopmental disorders. \
                   Journal of biochemistry, 163(6), 447–455. https://doi.org/10.1093/jb/mvy022")
    citation.write("**3.** Tippens, A. L., & Lee, A. (2007). Caldendrin, a neuron-specific modulator of Cav/1.2\
                    (L-type) Ca2+ channels. The Journal of biological chemistry, 282(11), 8464–8473. \
                   https://doi.org/10.1074/jbc.M611384200")
    
def plot_interpro(data = data, psd = psd, prob = prob, ID_count = ID_count):
    st.title('InterPro domains of EF-hand 1 and 2 binders')
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
    fig.update_layout(width=1200,height=800)

    st.plotly_chart(fig, theme='streamlit', use_container_width=True)

    binder_a_cont = st.container()
    binder_a_cont.write("Using BioGrid protein interaction database, we identified the interactors \
                        (binder B) of all EF1 and EF2 hand proteins (binder A) and identified all the \
                        potential motifs in their individual sequences. As a first step of this protein-protein \
                        interaction analysis, we were interested to see what is the general proportion of \
                        motif distribution among the EF1 and Ef2 hand Caldendrin-like proteins (binder A). \
                        The barplot represents the InterPro descriptions of all the identified motifs  from binder \
                        A sequences, plotted against the number of proteins consisting them. The data suggests that \
                        SRC3 domain, a well known protein interaction module, is present in around 45 binder A \
                        proteins making it the top motif observed in this dataset. Following SRC3, the PH like \
                        domain and PDZ domain are also present in similar proportions. The Pleckstrin homology (PH) \
                        like domains and PDZ domains are usually present in cytoskeleton to support the cellular architecture \
                        and also PH like domain is involved in intracellular signalling.")
    



########### Plot2 ###########
#Biomart query

    st.title('InterPro domains of post-synaptic density proteins')
    n_slider_2 = st.slider('Select the minimum value of n:', min_value=0, max_value=int(psd['n_PSD'].max()), value=100)
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
    fig2.update_layout(width=1200,height=800)
    st.plotly_chart(fig2, theme='streamlit', use_container_width=True)
    PSD_domain_cont = st.container()
    PSD_domain_cont.write("We used the PSD protein dataset to look \
                          at the distribution of each motif identified by InterPro. \
                          Like Figure 1, Figure 2 is a barplot of number of proteins \
                          consisting of each motif against the InterPro ID of each motif. \
                          The data suggests that PH-like domain and SRC3 domains are present \
                          in most of the PSD proteins which includes the EF1 and EF2 hand \
                          proteins (binder A proteins)")
    ########### Plot3 ###########
    #correlation plot

    st.title('Proportion of Interpro domains \ncomparing EF-hand 1 & 2 binders vs. total post synaptic density proteins\n')
    #calculate the proportion of each InterPro domain in the binder set
    
    search_term = st.text_input('Search for a data point from IPR number:')
    if search_term:
        prob = prob[prob.InterPro.str.contains(search_term, regex= True, na=False)]

    fig3 = px.scatter(prob, x = 'prop_x', 
                    y = 'prop_y', 
                    hover_data = ['Interpro.Description'], 
                    color = 'prop_x')

    fig3.update_layout(width=1200,height=800)
    fig3.update_xaxes(title_text='Proportion in EF-hand 1 & 2 binders', title_font={"size": 20}, title_standoff=25)
    fig3.update_yaxes(title_text='Proportion in PSD proteins', title_font={"size": 20}, title_standoff=25)
    st.plotly_chart(fig3, theme='streamlit', use_container_width=True)

    #######Interactors########
   
    #interactive slider
    n_slider_3 = st.slider('Select the minimum value of n:', min_value=0, max_value=int(ID_count['Count'].max()), value=50)
    ID_count_fil = ID_count.query("Count >= @n_slider_3")

    #interactive search
    search_term_2 = st.text_input('Search for domains based on name:')
    if search_term_2:
        ID_count_fil = ID_count[ID_count.Domains.str.contains(search_term_2, regex= True, na=False)]
    
    fig4 = px.bar(
        ID_count_fil,
        x='Domains',
        y='Count',
        title='Pfam Domains of Post-synaptic density proteins',
        hover_data=['Domains'],
        color="Count"
    )
    st.plotly_chart(fig4, theme='streamlit', use_container_width=True)

def disordered():
    st.header('Hidden Markov Model structure analysis')
    s = st.container()
    s.write("Displayed beneath is the most common schema of Caldendrin's disordered region identified by PHMMER")
    st.image('Disordered_schema.png', width = 800)
    st.divider()
    st.subheader('N-terminal Disordered region of Caldendrin frequently appears along double P-fam EF-hand 7 motif')
    text_cont1 = st.container()
    text_cont1.write("Consisting of ~ 200 amino acids, Caldendrin's disordered region \
            is a large component of the protein whose interaction preference is unknown. \
            A hidden markov model search for this architecture using EBI's HMMER suggestes that \
            this stretch of sequence is most commonly found in Caldendrin-like proteins across \
            the tree of life, with a strong presence in eukaryotes and birds")
    
    components.html('<iframe width="700" height="350" src="https://www.ebi.ac.uk/Tools/hmmer/results/DB852D26-D65D-11ED-A67D-A515B769042E.1/taxonomy/" title=“embedded website"></iframe>',
                     height=400, width=800)
    st.subheader('Opinion')
    comm_desc = st.container()
    comm_desc.write('Sites of significant functional protein-protein interactions are often highly conserved. \
                    We argue that this side of Caldendrin likely does not interacting with binding partners, since the MSA suggested fairly low \
                    sequence conservation compared to EF-hand motif. Nonetheless, the constant presence of this disordered stretch in Caldendrin-like proteins \
                    could still play other important functions not yet investigated.')

def MSA():
    st.header('Multiple Sequence Alignment')
    st.subheader('BLAST-ing against NCBI experimentally clustered data')
    msa_desc = st.container()
    msa_desc.write("Caldendrin is a ubiquitous protein found in many organisms, with fairly conserved N-terminal EF-hands cation binding motifs \
                   and a disordered N-terminal region.")
    msa_desc.write("To get a glimpse of the conservation of Caldendrin-like transcripts across the tree of life, \
                   we BLAST the human cannonical CABP1 sequence against NCBI's experimentally clustered database, \
                   then perform multiple sequence alignment using MUSCLE. \
                   Displayed beneath is a Dash-bio.AlignmentChart representation of the top few hundred hits, with sequence coverage ranging from 100% to 37% \
                   and E-values ranging from 6e-125 to 0.037")
    msa_desc.write('To get an overview of the alignment and modify the display range, scroll down to the bottom of the page adjust the slider. \
                   For a full-screen view, visit http://elpea.pythonanywhere.com')

    #https://blast.ncbi.nlm.nih.gov/Blast.cgi
    components.html('<iframe width="700" height="550" src="http://elpea.pythonanywhere.com/"></iframe>', height=800, width=800) 

    st.subheader('Motif from ligand B non-duplicated sequence')
    meme_cont = st.container()
    meme_cont.write("Running MEME motif sequence finder on the cannonical sequence of the ligand protein \
                    to search whether there is a sequence preference for proteins that bind EF-hand generally. \
                    Of all the sequences identified, the most common motif is displayed beneath, which occurs at 467 sequence \
                    out of the ~800 used for input.")
    st.image(image_motif, width = 800)

def known_interactors(interactor_df=interactor_df):
    st.header('Known Caldendrin binding partners')
    st.subheader('A. Ligand B proteins reported in the literature')

    fig5 = px.sunburst(
        interactor_df,
        path = ['Location', 'Protein', 'InterPro'],
        title='Caldendrin binding partners and associated domain',
        color='Location'
    )
    st.plotly_chart(fig5, theme='streamlit', use_container_width=True)
    
    partner_annotation = st.container()
    partner_annotation.write('The inner-most circle of the above chart represents the Caldendrin \
                             domains that these partners bind to, \
                             the middle circle represents the protein, \
                             and the outer-most circle are their associated InterPro ID. \
                             Not all listed IPR ID are involved in Caldendrin binding.')
    
    st.subheader('B. Comparing the reported domains with others that bind EF-hand motifs')
    
    color_plot =  px.scatter(prob, x = 'prop_x', 
                    y = 'prop_y', 
                    hover_data = ['Interpro.Description', 'InterPro', 'Binder_Protein_name'],
                    color = 'color_marker', 
                    labels={'color_marker': 'Classification', 
                            'prop_x': 'P(EF-hand binders) ', 
                            'prop_y': 'P(PSD protein)'},
                    color_discrete_sequence=px.colors.qualitative.Alphabet, 
                    color_discrete_map={'Other': '#0068C9', 'Known_binder': 'red'})
    #add dropdown menu to toggle log scale
    color_plot.update_traces(marker = dict(size =6,  line= dict(color = 'white', width = 1)))
    color_plot.add_shape(type="line", x0=0, y0=0, x1=0.04, y1=0.04, line=dict(color="grey", width=1, dash="dot"))
    

    st.plotly_chart(color_plot, theme='streamlit', use_container_width=True)
    
    plot_note = st.container()
    note = st.container()
    plot_note.write("The dash line is for x = y. The red colored scatter do not group in any particular cluster")
    
    note.write('The scatter plot shown above is a bit incomplete in describing the data density, \
               which is very strongly skewed towards smaller X and Y values.\
               Toggle the button beneath to see a 2d bin heatmap.')
    
    note_button = st.button('Extra heatmap')
    if note_button:
        test = px.density_heatmap(prob.query('prop_x < 0.0025').query('prop_y<0.0025'), 
                                  x = 'prop_x', y = 'prop_y', 
                                  color_continuous_scale='viridis', 
                                  labels = {'prop_x': 'P(EF-hand binders) ',
                                          'prop_y': 'P(PSD protein)'})
        st.plotly_chart(test, 
                        theme='streamlit',
                        use_container_width=True)
        
    #minor note

def discussion():
    st.header('Discussion')
    st.subheader('Initial thoughts')
    thoughts_cont = st.container()
    thoughts_cont.write("Our first attempt at this project is to create a set of **binder A** via structural search to identify\
                        proteins that are similar to Caldendrin. This approach is not viable since most human Caldendrin-like \
                        structures and sequences that we were able to find are just Caldendrin across different species, with very limited annotations on what \
                        their binding partners are. Thus, we decided to switch to motif querying. \
                        Since nature works in modules and motifs, we reason that looking at motifs which frequently interacts \
                        with EF-hands - Caldendrin's main structurally definitive component - can give us clues on \
                        potential interactors might be.")
    
    st.subheader('P(EF-hand binders) vs P(PSD protein)')
    thoughts_cont_2 = st.container()
    thoughts_cont_2.write(" We are treating *P(PSD-protein)*, the frequency at which given InterPro domain is found \
                          among PSD proteins, as a proxy for the commonality of the domain in general. \
                          *P(EF-hand binders)*, the frequency of the InterPro domain among EF-hand binders, is a variable \
                          that might reflect the prefrence of EF-hand for said InterPro domain. Since biologically relevant \
                          protein-protein interactions are often specific, we expect to see true EF-hand interactions  \
                          less frequently in the background (low **Y** frequency) and more frequently among positive hits for interacting partners (high **X** frequency). \
                          We reason that InterPro motifs that fall beneath the x = y as being more likely to represent candidate proteins. \
                          The probability of being a candidate partner would increase with")
    st.latex(r'''\left| X - Y \right|''')
    
    subthoughts_cont = st.container()
    subthoughts_cont.write("Most of Caldendrin's known binding candidates do fall under the x = y line, \
                           although the same can be said for the majority of our dataset. This is overall not too surprising, \
                           since EF-hand binders forms a subset of the collective PSD protein. The outlier of this trend is Cortactin, \
                           whose InterPro annotations are both common among the background proteins and EF-hand binders.")

    thoughts_cont_3 = st.container()
    thoughts_cont_3.write("Our heuristic to finding binders have several drawbacks. First, we are not taking into consideration the \
                        presence of orthogonal mutations. As proteins evolve and diverge, they can acquire mutations that \
                        enable them to bind to new ligands and no longer interact with the original ligand while still maintaining \
                        their motif. In addition, we do not have a way to asssess the confidence of our predictions statistically. \
                        Last but not least, we have not validate whether our heuristic is viable using available experimental data.")

    st.subheader('Future directions')
    future_cont = st.container()
    future_cont.write("A good improvement to the model is to also take into consideration the number of proteins that carries the give motif, both at the \
                      PSD and proteome level. More careful vetting of the data, as well as *ab initio* docking calculations and MD simulations of predicted targets\
                      interaction with Caldendrin could be useful in improving our predictions.")
    
    st.subheader('Closing remarks')
    st.markdown("Thank you professor Marcotte for giving us a comprehensive introduction to the wonderful world of Bioinformatics and Systems biology, for always keeping \
                      lectures interesting with the latest scientific update, and for inspiring us be curious scientists. \
                      We would also like to thank our TA Matt, for providing great coding advices, for being super helpful in explaining concepts, and for introducing us to \
                      Streamlit, the backbone of this website.")
    st.markdown("Code and data can be found at https://github.com/Lisa-Phan/Bioinfo_Project_Spring2023/tree/Lisa-Phan-app-v1")

    
    
    



# Define dropdown menu options
options = ['Introduction', 'MSA and motif', 
           'InterPro domains distribution',
           'Disordered domain', 'Known interactors', 
           'Discussion']

# Create dropdown menu
choice = st.sidebar.selectbox('Choose an option:', options)

# Display selected option
if choice == 'Introduction':
    intro()
elif choice == 'InterPro domains distribution':
    plot_interpro()
elif choice == 'Disordered domain':
    disordered()
elif choice == 'MSA and motif':
    MSA()
elif choice == 'Known interactors':
    known_interactors()
elif choice == 'Discussion':
    discussion()






