library(shiny)
library(tidyverse)
library(viridis)
library(shinydashboard)
library(dashboardthemes)


server <- function(input, output) {
  
  
  
  output$text1 <- renderText({
    "Difference in gene set mutation rate versus expected mutation rate. Substitution mutations determined from whole genome sequencing of human induced pluripotent stem cells (HiPSC) treated with environmental chemicals."
  })
  
  # no longer using text2
  output$text2 <- renderText({
    "Per gene mutation rates of selected disease gene sets. Substitution mutations determined from whole genome sequencing of human induced pluripotent stem cells (HiPSC) treated with environmental chemicals. Grey points indicate mutation rates in control treatment conditions. Red lines indicate expected mutation rates for each chemical treatment, calculated with a bootstraping method."
  })
  
  output$text3 <- renderText({
    "Input gene names filtered to keep only those included in GENCODE release 37."
  })
  
  
  
# filter the user input genes to keep only genes with valid grch37 gene names

  user_grch37_genes <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    req(input$file1, file.exists(input$file1$datapath))
    req(input$select)
    input_genes = read_csv(input$file1$datapath) %>% janitor::clean_names()

    grch37_gene_names = read_csv("./grch37_gene_names.csv") %>% pull(gene_name)
    
    user_genes = lapply(input_genes,  function(x) {x[x %in% grch37_gene_names]})
    user_genes
    
    # in the list of lists above, each list is not the same length, so we cannot convert to data frame for rendering. make them all the same length by adding NAs to the end of each list
    #get the length of the longest list
    maxlen <- max(lengths(user_genes))
    #add NAs to make each list as long as the longest list
    user_genes2 <- lapply(user_genes, function(lst) c(lst, rep(NA, maxlen - length(lst))))
    
    as.data.frame(user_genes2)
 
  })
  
  
# Gene body data for heatmap    
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
    grch37_gene_names = read_csv("./grch37_gene_names.csv")
    
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
        
        gene_number = gene_list_x %>% filter(gene %in% grch37_gene_names$gene_name)
        
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
                    genes = length(gene_number$gene),
                    wells = length(n)) %>% 
          mutate(gene_wells = genes*wells) 
        
        return(data1)}
      input_mutations = input_mutations %>% bind_cols(chem_list_df) %>% 
        mutate(name = name) 
    }
    
    input = input_list %>% 
      dplyr::rename(input_mutations = total_mutations,input_genewells = gene_wells, chemical = chem_list) %>% 
      filter(chemical != "Control")%>% 
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
          
          test1 = binom.test(gene_list_x[["input_mutations"]][[row]],
                             gene_list_x[["input_genewells"]][[row]],
                             gene_list_x[["mean"]][[row]],
                             alternative = "two.sided")
          broom::tidy(test1) %>% mutate(chemical = gene_list_x[["chemical"]][[row]])
          
        }
        #happens between the rbinds
        input_tests = input_tests %>% mutate(name = list_name)
      }
    
    
    
    # Add back bin_test_data to add back the control rate (mean) and calucalte differences and fold changes
    input_binomial_tests_final = input_binomial_tests %>% 
      mutate(treatment = "Treated") %>% 
      left_join(bin_test_data) %>% 
      mutate(difference = estimate - mean,
             fold_change = estimate/mean)
    
    
    #arrange in order of highest to lowest estimate, based on radiation
    input_binomial_tests_final$chemical <- factor(input_binomial_tests_final$chemical,levels = c("Radiation", "PAH", "Alkylating Agent", "Nitro-PAH", "ROS/NOS", "Heterocyclic Amine", "Drug Therapy", "DNA Damage Response Inhibitors", "Aromatic Amine", "Metal", "Nitrosamine", "Aldehydes", "Other"))
    
    
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
    grch37_gene_names = read_csv("./grch37_gene_names.csv")
    
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
        
        gene_number = gene_list_x %>% filter(gene %in% grch37_gene_names$gene_name)
        
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
                    genes = length(gene_number$gene),
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
      filter(chemical != "Control")%>% 
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
          
          test1 = binom.test(gene_list_x[["input_mutations"]][[row]],
                             gene_list_x[["input_genewells"]][[row]],
                             gene_list_x[["mean"]][[row]],
                             alternative = "two.sided")
          broom::tidy(test1) %>% mutate(chemical = gene_list_x[["chemical"]][[row]])
          
        }
        #happens between the rbinds
        input_tests = input_tests %>% mutate(name = list_name)
      }
    
    
    
    input_binomial_tests_final = input_binomial_tests %>% 
      mutate(treatment = "Treated") %>% 
      left_join(bin_test_data) %>% 
      mutate(difference = estimate - mean,
             fold_change = estimate/mean)
    
    input_binomial_tests_final$chemical <- factor(input_binomial_tests_final$chemical,levels = c("Radiation", "PAH", "Alkylating Agent", "Nitro-PAH", "ROS/NOS", "Heterocyclic Amine", "Drug Therapy", "DNA Damage Response Inhibitors", "Aromatic Amine", "Metal", "Nitrosamine", "Aldehydes", "Other"))
    
    
    
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
      write_csv(heatmap_data(), file)
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
      mutate(significance = ifelse(adj_p > 0.05, significance, "**"))%>% 
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
  
  
  heatmap_data <- reactive({
    
    final_data = final_data() %>% 
      mutate(adj_p = p.adjust(p.value, method = "bonferroni")) %>%
      mutate(significance = ifelse(p.value > 0.05, NA, "*"))%>%
      mutate(significance = ifelse(adj_p > 0.05, significance, "**"))%>% 
      mutate(significance = ifelse(significance =="*" & difference<0, "???",significance))%>% 
      mutate(significance = ifelse(significance =="**" & difference<0, "??????",significance))%>%
      mutate(color = ifelse(difference>0 & !is.na(significance),"darkturquoise", NA ))%>% 
      mutate(color = ifelse(difference<0 & !is.na(significance),"red", color))%>% 
      mutate(color = ifelse(is.na(significance),"black", color))
    
   final_data
    
  })
  
#the example input in the left panel  
  output$example <- renderTable({
    example = read_csv("./test_input.csv")
    example})    
#the users gene lists are filtered on grch37 gene names. The filtered genes are shown back to the user  
  output$contents3 <- renderTable({
    req(user_grch37_genes())
    user_grch37_genes()})

#data in the results table tab    
    output$contents2 <- renderTable({
    req(heatmap_data())
    heatmap_data()})
#old plot not shown anymore  
  output$contents <- renderPlot({
    req(plot1())
    plot1()})
#heatmap plot  
  output$heatmap_gene <- renderPlot({
    req(heatmap_gene())
    heatmap_gene()
    
    
    
    
    
  })
}
