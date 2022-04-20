library(shiny)
library(tidyverse)
library(viridis)
library(shinydashboard)
library(dashboardthemes)


server <- function(input, output) {
  
  
  
  output$text1 <- renderText({
    "Fold-change in gene set mutation rate over expected mutation rate. Substitution mutations determined from whole genome sequencing of human induced pluripotent stem cells (HiPSC) treated with environmental chemicals."
  })
  
  
  output$text2 <- renderText({
    "Per gene mutation rates of selected disease gene sets. Substitution mutations determined from whole genome sequencing of human induced pluripotent stem cells (HiPSC) treated with environmental chemicals. Grey points indicate mutation rates in control treatment conditions. Red lines indicate expected mutation rates for each chemical treatment, calculated with a bootstraping method."
  })
  
  
  final_data_gene_body <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    req(input$file1, file.exists(input$file1$datapath))
    req(input$select)
    input_genes = read_csv(input$file1$datapath)
    
    library(tidyverse)
    library(foreach)
    
    control_props = read_csv("./control_props.csv")
    
    wells = read_csv("./wells.csv")
    chem_list = wells %>% ungroup() %>%  distinct(group) %>% pull(group)
    subs = read_csv("./subs_gene_body.csv")
    hline_data = read_csv("./hline_data.csv")
    all_gene_length = read_csv("./all_gene_length.csv")
    
    #loop
    #input_genes=input_genes %>% dplyr::select(CHD) %>% #dplyr::rename(gene=CHD)
    
    gene_list_names = colnames(input_genes)
    #####################
    ####################
    ###### setup loop for binomial tests
    #####################
    ######################
    input_list = foreach(name = gene_list_names, .combine = "rbind") %do% {
      gene_list_x = input_genes %>% dplyr::select(name) %>% dplyr::rename(gene=name) %>% filter(!is.na(gene))
      
      
      total_nucleotides = gene_list_x %>% left_join(all_gene_length) %>% summarize(mean1 = mean(Length, na.rm=TRUE)) %>% pull(mean1)
      
      chem_list_df = data.frame(chem_list)
      input_mutations = foreach(chemical = chem_list, .combine = "rbind") %do% {
        
        data1 = gene_list_x  %>% 
          left_join(subs) %>%    
          filter(group == chemical) %>% 
          group_by(group,chem,dose,replicate) %>% 
          summarize(n=n()) %>% 
          full_join(wells) %>% 
          filter(group == chemical) %>% 
          mutate(n = replace_na(n, 0)) %>% 
          ungroup() %>% 
          summarize(total_mutations = sum(n),
                    genes = length(gene_list_x$gene),
                    wells = length(n)) %>% 
          mutate(gene_wells = genes*wells) %>% 
          mutate(gene_wells_nucleotides = gene_wells*total_nucleotides)
        
        return(data1)}
      input_mutations = input_mutations %>% bind_cols(chem_list_df) %>% 
        mutate(name = name) 
    }
    
    #chem_list_df = data.frame(chem_list)
    #input_mutations$chemical = chem_list_df$chem_list
    
    input = input_list %>% 
      dplyr::rename(input_mutations = total_mutations,input_genewells = gene_wells, chemical = chem_list) %>% 
      dplyr::select(chemical, input_mutations, input_genewells, name)
    
    bin_test_data = input %>% left_join(control_props)
    
    ##################
    ###################
    #Binomial test loop
    ###################
    #################
    input_binomial_tests = foreach (list_name = gene_list_names, .combine = "rbind") %do%
      {
        gene_list_x = bin_test_data %>% filter(name == list_name)
        
        input_tests=foreach (row = 1:nrow(gene_list_x), .combine = "rbind") %do% {
          
          test1 = binom.test(gene_list_x[["input_mutations"]][[row]],gene_list_x[["input_genewells"]][[row]],gene_list_x[["mean"]][[row]])
          broom::tidy(test1) %>% mutate(chemical = gene_list_x[["chemical"]][[row]])
          
        }
        #happens between the rbinds
        input_tests = input_tests %>% mutate(name = list_name)
      }
    
    
    
    ################
    ################
    #add control loop
    ##################
    ################
    input_binomial_tests2 = input_binomial_tests %>% mutate(estimate2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$estimate2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$estimate
      
    }
    
    #confhigh
    input_binomial_tests2 = input_binomial_tests2 %>% mutate(conf.high2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$conf.high2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$conf.high
      
    }
    
    #conflow
    input_binomial_tests2 = input_binomial_tests2 %>% mutate(conf.low2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$conf.low2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$conf.low
      
    }
    
    
    input_binomial_tests2 = input_binomial_tests2 %>% 
      dplyr::select(-estimate, -conf.high, -conf.low) %>% 
      dplyr::rename(estimate=estimate2,
                    conf.high = conf.high2,
                    conf.low = conf.low2) %>% 
      mutate(treatment = "Control")
    
    input_binomial_tests = input_binomial_tests %>% 
      mutate(treatment = "Treated")
    
    input_binomial_tests_final = input_binomial_tests %>% 
      bind_rows(input_binomial_tests2) %>% 
      filter(chemical != "Control")
    
    ############
    ###############
    ##### add hline 
    ###############
    ##############
    ############hlines for expected mutations for random gene sets
    Radiation_hline = control_props %>% filter(chemical == "Radiation") %>% pull(mean)
    PAH_hline = control_props %>% filter(chemical == "PAH")%>% pull(mean)
    Alkylating_Agent_hline = control_props %>% filter(chemical == "Alkylating Agent")%>% pull(mean)
    NitroPAH_hline = control_props %>% filter(chemical == "Nitro-PAH")%>% pull(mean)
    ROSNOS_hline = control_props %>% filter(chemical == "ROS/NOS")%>% pull(mean)
    Heterocyclic_Amine_hline = control_props %>% filter(chemical == "Heterocyclic Amine")%>% pull(mean)
    Drug_Therapy_hline = control_props %>% filter(chemical == "Drug Therapy")%>% pull(mean)
    DNA_Damage_Response_Inhibitors_hline = control_props %>% filter(chemical == "DNA Damage Response Inhibitors")%>% pull(mean)
    Aromatic_Amine_hline = control_props %>% filter(chemical == "Aromatic Amine")%>% pull(mean)
    Metal_hline = control_props %>% filter(chemical == "Metal")%>% pull(mean)
    Nitrosamine_hline = control_props %>% filter(chemical == "Nitrosamine")%>% pull(mean)
    Aldehydes_hline = control_props %>% filter(chemical == "Aldehydes")%>% pull(mean)
    Other_hline = control_props %>% filter(chemical == "Other")%>% pull(mean)
    Control_hline = control_props %>% filter(chemical == "Control")%>% pull(mean)
    
    
    
    
    #write_csv(hline_data, "./hline_data.csv")
    ##############add fold change
    
    #final_data$gene_set <- factor(final_data$gene_set,levels = c("Input genes" ,"Schizophrenia", "Autism", "Obesity", "ADHD", "Coronary artery disease", "Alzheimers", "ALS", "CHD", "T2 Diabetes", "Oral cleft"))
    
    
    
    input_binomial_tests_final = input_binomial_tests_final %>% mutate(fold_change = case_when(chemical == "Aldehydes" ~ estimate/Aldehydes_hline,
                                                                                               chemical == "Alkylating Agent" ~ estimate/Alkylating_Agent_hline,
                                                                                               chemical == "Aromatic Amine" ~ estimate/Aromatic_Amine_hline,
                                                                                               chemical == "Control" ~ estimate/Control_hline,
                                                                                               chemical == "DNA Damage Response Inhibitors" ~ estimate/DNA_Damage_Response_Inhibitors_hline,
                                                                                               chemical == "Heterocyclic Amine" ~ estimate/Heterocyclic_Amine_hline,
                                                                                               chemical == "Metal" ~ estimate/Metal_hline,
                                                                                               chemical == "Nitro-PAH" ~ estimate/NitroPAH_hline,
                                                                                               chemical == "Nitrosamine" ~ estimate/Nitrosamine_hline,
                                                                                               chemical == "Other" ~ estimate/Other_hline,
                                                                                               chemical == "PAH" ~ estimate/PAH_hline,
                                                                                               chemical == "Radiation" ~ estimate/Radiation_hline,
                                                                                               chemical == "ROS/NOS" ~ estimate/ROSNOS_hline,
                                                                                               chemical == "Drug Therapy" ~ estimate/Drug_Therapy_hline))
    
    input_binomial_tests_final = input_binomial_tests_final %>% mutate(difference = case_when(chemical == "Aldehydes" ~ estimate-Aldehydes_hline,
                                                                                              chemical == "Alkylating Agent" ~ estimate-Alkylating_Agent_hline,
                                                                                              chemical == "Aromatic Amine" ~ estimate-Aromatic_Amine_hline,
                                                                                              chemical == "Control" ~ estimate-Control_hline,
                                                                                              chemical == "DNA Damage Response Inhibitors" ~ estimate-DNA_Damage_Response_Inhibitors_hline,
                                                                                              chemical == "Heterocyclic Amine" ~ estimate-Heterocyclic_Amine_hline,
                                                                                              chemical == "Metal" ~ estimate-Metal_hline,
                                                                                              chemical == "Nitro-PAH" ~ estimate-NitroPAH_hline,
                                                                                              chemical == "Nitrosamine" ~ estimate-Nitrosamine_hline,
                                                                                              chemical == "Other" ~ estimate-Other_hline,
                                                                                              chemical == "PAH" ~ estimate-PAH_hline,
                                                                                              chemical == "Radiation" ~ estimate-Radiation_hline,
                                                                                              chemical == "ROS/NOS" ~ estimate-ROSNOS_hline,
                                                                                              chemical == "Drug Therapy" ~ estimate-Drug_Therapy_hline))
    
    #arrange in order of highest to lowest estimate, based on radiation
    order = input_binomial_tests_final %>% filter(chemical == "Radiation",
                                                  treatment=="Treated") %>% 
      arrange(-estimate) %>% pull(name)
    #re order factors based on above
    input_binomial_tests_final$name <- factor(input_binomial_tests_final$name,levels = order)
    
    input_binomial_tests_final$chemical <- factor(input_binomial_tests_final$chemical,levels = c("Radiation", "PAH", "Alkylating Agent", "Nitro-PAH", "ROS/NOS", "Heterocyclic Amine", "Drug Therapy", "DNA Damage Response Inhibitors", "Aromatic Amine", "Metal", "Nitrosamine", "Aldehydes", "Other", "Control"))
    
    
    
    
    input_binomial_tests_final
    
  })
  ############
  ################
  #################
  ###############
  #################
  #################
  ################
  ################
  ###############
  
  final_data_cds = reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    req(input$file1, file.exists(input$file1$datapath))
    req(input$select)
    input_genes = read_csv(input$file1$datapath)
    
    library(tidyverse)
    library(foreach)
    
    control_props = read_csv("./control_props_cds.csv")
    wells = read_csv("./wells.csv")
    chem_list = wells %>% ungroup() %>%  distinct(group) %>% pull(group)
    subs = read_csv("./subs_cds.csv")
    hline_data = read_csv("./hline_data_cds.csv")
    all_gene_length = read_csv("./all_gene_length_coding.csv")
    
    #loop
    #input_genes=input_genes %>% dplyr::select(CHD) %>% #dplyr::rename(gene=CHD)
    
    gene_list_names = colnames(input_genes)
    #####################
    ####################
    ###### setup loop for binomial tests
    #####################
    ######################
    input_list = foreach(name = gene_list_names, .combine = "rbind") %do% {
      gene_list_x = input_genes %>% dplyr::select(name) %>% dplyr::rename(gene=name) %>% filter(!is.na(gene))
      
      chem_list_df = data.frame(chem_list)
      input_mutations = foreach(chemical = chem_list, .combine = "rbind") %do% {
        
        data1 = gene_list_x  %>% 
          left_join(subs) %>%    
          filter(group == chemical) %>% 
          group_by(group,chem,dose,replicate) %>% 
          summarize(n=n()) %>% 
          full_join(wells) %>% 
          filter(group == chemical) %>% 
          mutate(n = replace_na(n, 0)) %>% 
          ungroup() %>% 
          summarize(total_mutations = sum(n),
                    genes = length(gene_list_x$gene),
                    wells = length(n)) %>% 
          mutate(gene_wells = genes*wells) 
        
        return(data1)}
      input_mutations = input_mutations %>% bind_cols(chem_list_df) %>% 
        mutate(name = name) 
    }
    
    #chem_list_df = data.frame(chem_list)
    #input_mutations$chemical = chem_list_df$chem_list
    
    input = input_list %>% 
      dplyr::rename(input_mutations = total_mutations,input_genewells = gene_wells, chemical = chem_list) %>% 
      dplyr::select(chemical, input_mutations, input_genewells, name)
    
    bin_test_data = input %>% left_join(control_props)
    
    ##################
    ###################
    #Binomial test loop
    ###################
    #################
    input_binomial_tests = foreach (list_name = gene_list_names, .combine = "rbind") %do%
      {
        gene_list_x = bin_test_data %>% filter(name == list_name)
        
        input_tests=foreach (row = 1:nrow(gene_list_x), .combine = "rbind") %do% {
          
          test1 = binom.test(gene_list_x[["input_mutations"]][[row]],gene_list_x[["input_genewells"]][[row]],gene_list_x[["mean"]][[row]])
          broom::tidy(test1) %>% mutate(chemical = gene_list_x[["chemical"]][[row]])
          
        }
        #happens between the rbinds
        input_tests = input_tests %>% mutate(name = list_name)
      }
    
    
    
    ################
    ################
    #add control loop
    ##################
    ################
    input_binomial_tests2 = input_binomial_tests %>% mutate(estimate2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$estimate2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$estimate
      
    }
    
    #confhigh
    input_binomial_tests2 = input_binomial_tests2 %>% mutate(conf.high2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$conf.high2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$conf.high
      
    }
    
    #conflow
    input_binomial_tests2 = input_binomial_tests2 %>% mutate(conf.low2 = 0)
    for (i in gene_list_names) {
      input_binomial_tests2[input_binomial_tests2$name == i,]$conf.low2 <- input_binomial_tests2[input_binomial_tests2$chemical == "Control"&input_binomial_tests2$name == i,]$conf.low
      
    }
    
    
    input_binomial_tests2 = input_binomial_tests2 %>% 
      dplyr::select(-estimate, -conf.high, -conf.low) %>% 
      dplyr::rename(estimate=estimate2,
                    conf.high = conf.high2,
                    conf.low = conf.low2) %>% 
      mutate(treatment = "Control")
    
    input_binomial_tests = input_binomial_tests %>% 
      mutate(treatment = "Treated")
    
    input_binomial_tests_final = input_binomial_tests %>% 
      bind_rows(input_binomial_tests2) %>% 
      filter(chemical != "Control")
    
    ############
    ###############
    ##### add hline 
    ###############
    ##############
    ############hlines for expected mutations for random gene sets
    Radiation_hline = control_props %>% filter(chemical == "Radiation") %>% pull(mean)
    PAH_hline = control_props %>% filter(chemical == "PAH")%>% pull(mean)
    Alkylating_Agent_hline = control_props %>% filter(chemical == "Alkylating Agent")%>% pull(mean)
    NitroPAH_hline = control_props %>% filter(chemical == "Nitro-PAH")%>% pull(mean)
    ROSNOS_hline = control_props %>% filter(chemical == "ROS/NOS")%>% pull(mean)
    Heterocyclic_Amine_hline = control_props %>% filter(chemical == "Heterocyclic Amine")%>% pull(mean)
    Drug_Therapy_hline = control_props %>% filter(chemical == "Drug Therapy")%>% pull(mean)
    DNA_Damage_Response_Inhibitors_hline = control_props %>% filter(chemical == "DNA Damage Response Inhibitors")%>% pull(mean)
    Aromatic_Amine_hline = control_props %>% filter(chemical == "Aromatic Amine")%>% pull(mean)
    Metal_hline = control_props %>% filter(chemical == "Metal")%>% pull(mean)
    Nitrosamine_hline = control_props %>% filter(chemical == "Nitrosamine")%>% pull(mean)
    Aldehydes_hline = control_props %>% filter(chemical == "Aldehydes")%>% pull(mean)
    Other_hline = control_props %>% filter(chemical == "Other")%>% pull(mean)
    Control_hline = control_props %>% filter(chemical == "Control")%>% pull(mean)
    
    
    
    
    #write_csv(hline_data, "./hline_data.csv")
    ##############add fold change
    
    #final_data$gene_set <- factor(final_data$gene_set,levels = c("Input genes" ,"Schizophrenia", "Autism", "Obesity", "ADHD", "Coronary artery disease", "Alzheimers", "ALS", "CHD", "T2 Diabetes", "Oral cleft"))
    
    
    
    input_binomial_tests_final = input_binomial_tests_final %>% mutate(fold_change = case_when(chemical == "Aldehydes" ~ estimate/Aldehydes_hline,
                                                                                               chemical == "Alkylating Agent" ~ estimate/Alkylating_Agent_hline,
                                                                                               chemical == "Aromatic Amine" ~ estimate/Aromatic_Amine_hline,
                                                                                               chemical == "Control" ~ estimate/Control_hline,
                                                                                               chemical == "DNA Damage Response Inhibitors" ~ estimate/DNA_Damage_Response_Inhibitors_hline,
                                                                                               chemical == "Heterocyclic Amine" ~ estimate/Heterocyclic_Amine_hline,
                                                                                               chemical == "Metal" ~ estimate/Metal_hline,
                                                                                               chemical == "Nitro-PAH" ~ estimate/NitroPAH_hline,
                                                                                               chemical == "Nitrosamine" ~ estimate/Nitrosamine_hline,
                                                                                               chemical == "Other" ~ estimate/Other_hline,
                                                                                               chemical == "PAH" ~ estimate/PAH_hline,
                                                                                               chemical == "Radiation" ~ estimate/Radiation_hline,
                                                                                               chemical == "ROS/NOS" ~ estimate/ROSNOS_hline,
                                                                                               chemical == "Drug Therapy" ~ estimate/Drug_Therapy_hline))
    
    
    
    input_binomial_tests_final = input_binomial_tests_final %>% mutate(difference = case_when(chemical == "Aldehydes" ~ estimate-Aldehydes_hline,
                                                                                              chemical == "Alkylating Agent" ~ estimate-Alkylating_Agent_hline,
                                                                                              chemical == "Aromatic Amine" ~ estimate-Aromatic_Amine_hline,
                                                                                              chemical == "Control" ~ estimate-Control_hline,
                                                                                              chemical == "DNA Damage Response Inhibitors" ~ estimate-DNA_Damage_Response_Inhibitors_hline,
                                                                                              chemical == "Heterocyclic Amine" ~ estimate-Heterocyclic_Amine_hline,
                                                                                              chemical == "Metal" ~ estimate-Metal_hline,
                                                                                              chemical == "Nitro-PAH" ~ estimate-NitroPAH_hline,
                                                                                              chemical == "Nitrosamine" ~ estimate-Nitrosamine_hline,
                                                                                              chemical == "Other" ~ estimate-Other_hline,
                                                                                              chemical == "PAH" ~ estimate-PAH_hline,
                                                                                              chemical == "Radiation" ~ estimate-Radiation_hline,
                                                                                              chemical == "ROS/NOS" ~ estimate-ROSNOS_hline,
                                                                                              chemical == "Drug Therapy" ~ estimate-Drug_Therapy_hline))
    
    #arrange in order of highest to lowest estimate, based on radiation
    order = input_binomial_tests_final %>% filter(chemical == "Radiation",
                                                  treatment=="Treated") %>% 
      arrange(-estimate) %>% pull(name)
    #re order factors based on above
    input_binomial_tests_final$name <- factor(input_binomial_tests_final$name,levels = order)
    
    input_binomial_tests_final$chemical <- factor(input_binomial_tests_final$chemical,levels = c("Radiation", "PAH", "Alkylating Agent", "Nitro-PAH", "ROS/NOS", "Heterocyclic Amine", "Drug Therapy", "DNA Damage Response Inhibitors", "Aromatic Amine", "Metal", "Nitrosamine", "Aldehydes", "Other", "Control"))
    
    #####cds
    
    
    input_binomial_tests_final
  })
  
  final_data = reactive({ 
    req(final_data_gene_body())
    if (input$select == '1') {
      data <- final_data_gene_body()
    }
    else
    {
      data <- final_data_cds()
    }
  })
  
  output$download1 <- downloadHandler(
    filename = function() {
      "results.csv"
    },
    content = function(file) {
      write_csv(final_data(), file)
    }
  )  
  
  plot1 <- reactive({
    
    if (input$select == '1') {
      hline_data = read_csv("./hline_data.csv")
    }
    else
    {
      hline_data = read_csv("./hline_data_cds.csv")
    }
    
    final_data = final_data()
    
    bbtheme <- theme(axis.text.x = element_text(size=15,face = "bold"),
                     axis.text.y = element_text(size=15, face = "bold"),
                     axis.title.x = element_text(size=20, face = "bold"),
                     axis.title.y = element_text(size=20, face = "bold"),
                     plot.title = element_text(size=20, face = "bold"))
    
    ggplot(data = final_data, aes(x = name, y = estimate, color = treatment))+
      geom_point(position=position_dodge(.2))+
      geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2,
                    position=position_dodge(.2))+
      geom_hline(data=hline_data, aes(yintercept = hline), color = "red")+
      labs(y = "Mutation rate",
           x = "Input gene set")+
      #geom_hline(yintercept = Aldehydes_hline, color = "red")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.y = element_blank())+
      scale_color_manual(values = c("gray50","black"))+
      facet_wrap(~factor(chemical, levels = c("Radiation", "PAH", "Alkylating Agent", "Nitro-PAH", "ROS/NOS", "Heterocyclic Amine", "Drug Therapy", "DNA Damage Response Inhibitors", "Aromatic Amine", "Metal", "Nitrosamine", "Aldehydes", "Other", "Control")))+
      bbtheme
  })
  
  
  ##################
  ###heatmap
  
  
  
  heatmap_gene <- reactive({
    
    final_data = final_data() %>% filter(treatment == "Treated") %>% 
      mutate(adj_p = p.adjust(p.value, method = "bonferroni")) %>%
      mutate(significance = ifelse(p.value > 0.05, NA, "*"))%>%
      mutate(significance = ifelse(adj_p > 0.05, significance, "***"))%>% 
      mutate(significance = ifelse(significance == "hi", NA, significance))
    
    bbtheme <- theme(axis.text.x = element_text(size=15,face = "bold"),
                     axis.text.y = element_text(size=15, face = "bold"),
                     axis.title.x = element_text(size=20, face = "bold"),
                     axis.title.y = element_text(size=20, face = "bold"),
                     plot.title = element_text(size=20, face = "bold"))
    
    ggplot(final_data, aes(x = name, y = chemical, fill = difference)) + 
      geom_tile()+
      geom_text(aes(label = significance), color = "darkturquoise", fontface = "bold", size = 16) +
      scale_fill_viridis(option = "B")+
      labs(fill = "Observed\nminus expected\nmutation rate",
           y = "Chemical class",
           x = "Input gene set")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      bbtheme
    
  })
  
  
  output$example <- renderTable({
    example = read_csv("./test_input.csv")
    example})    
  
  output$contents2 <- renderTable({
    req(final_data())
    final_data()})
  
  output$contents <- renderPlot({
    req(plot1())
    plot1()})
  
  output$heatmap_gene <- renderPlot({
    req(heatmap_gene())
    heatmap_gene()
    
    
    
    
    
  })
}
