##Attempting to make a really cool figure in R, the life cycle to replace
#the ugly one made from objects in power point


##Using pipes 
  
library(DiagrammeR)
library(magrittr)

#Example from DiagrammeR github page 
graph <-
  create_graph() %>%
  set_graph_name("DAG") %>%
  set_global_graph_attr("graph", "overlap", "true") %>%
  set_global_graph_attr("graph", "fixedsize", "true") %>%
  set_global_graph_attr("node", "color", "blue") %>%
  set_global_graph_attr("node", "fontname", "Helvetica") %>%
  add_n_nodes(11) %>%
  select_nodes_by_id(1:4) %>% 
  set_node_attr_with_selection("shape", "box") %>%
  set_node_attr_with_selection("type", "box") %>%
  clear_selection %>%
  select_nodes_by_id(5:7) %>% 
  set_node_attr_with_selection("shape", "circle") %>%
  set_node_attr_with_selection("type", "circle") %>%
  clear_selection %>%
  select_nodes_by_id(8:11) %>% 
  set_node_attr_with_selection("shape", "box") %>%
  set_node_attr_with_selection("type", "box") %>%
  clear_selection %>%
  add_edge(1, 5) %>% add_edge(2, 6) %>%
  add_edge(3, 9) %>% add_edge(4, 7) %>%
  add_edge(5, 8) %>% add_edge(5, 10) %>%
  add_edge(7, 11) %>% 
  select_edges %>%
  set_edge_attr_with_selection("color", "green") %>%
  add_edge(1, 8) %>% add_edge(3, 6) %>%
  add_edge(3, 11) %>% add_edge(3, 7) %>%
  add_edge(5, 9) %>% add_edge(6, 10) %>%
  select_edges("color", "^$") %>%
  set_edge_attr_with_selection("color", "red") %>%
  clear_selection

render_graph(graph)  

#Take a crack at it 
graph <-
  create_graph() %>%
  set_graph_name("DAG") %>%
  set_global_graph_attr("graph", "overlap", "true") %>%
  set_global_graph_attr("graph", "fixedsize", "true") %>%
  set_global_graph_attr("node", "color", "blue") %>%
  set_global_graph_attr("node", "fontname", "Helvetica") %>%
  add_n_nodes(4) %>%
  select_nodes_by_id(1:2) %>% 
  set_node_attr_with_selection("shape", "circle") %>%
  set_node_attr_with_selection("type", "circle") %>%
  clear_selection %>%
  select_nodes_by_id(3:4) %>% 
  set_node_attr_with_selection("shape", "circle") %>%
  set_node_attr_with_selection("type", "circle") %>%
  clear_selection %>%
  add_edge(1, 2) %>% add_edge(2, 1) %>%
  add_edge(3, 4) %>% add_edge(4, 3) %>%
  select_edges %>%
  set_edge_attr_with_selection("color", "green") %>%
  add_edge(1, 1) %>% add_edge(2, 2) %>%
  add_edge(3, 3) %>% add_edge(4, 4) %>%
  select_edges("color", "^$") %>%
  set_edge_attr_with_selection("color", "red") %>%
  clear_selection

render_graph(graph)  

####### Using a different package 
#Diagram 

# Create population matrix
#
Numgenerations <- 6
DiffMat <- matrix(data = 0, nrow = Numgenerations, ncol = Numgenerations)
AA <- as.data.frame(DiffMat)
AA[[1,4]] <- "f[3]"
AA[[1,5]] <- "f[4]"
AA[[1,6]] <- "f[5]"
#
AA[[2,1]] <- "s[list(0,1)]"
AA[[3,2]] <- "s[list(1,2)]"
AA[[4,3]] <- "s[list(2,3)]"
AA[[5,4]] <- "s[list(3,4)]"
AA[[6,5]] <- "s[list(4,5)]"
#
name <- c(expression(Age[0]), expression(Age[1]), expression(Age[2]),
               expression(Age[3]), expression(Age[4]), expression(Age[5]))

plotmat(A = AA, pos = 6, curve = 0.7, name = name, lwd = 2,
             arr.len = 0.6, arr.width = 0.25, my = -0.2,
             box.size = 0.05, arr.type = "triangle", dtext = 0.95,
             main = "Age-structured population model 1")
