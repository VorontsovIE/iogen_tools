require_relative 'cluster_formatter'

class ClusterNewickFormatter < ClusterFormatter
  def create_newick_html(branch_len_meth, size={})
    newick_inner_content = create_newick_inner_content(branch_len_meth)
    size_x = size[:x] || 2000
    size_y = size[:y] || 2000
    
    <<-RESULT
    <html>
    <head>
      <script type="text/javascript" src="js/raphael/raphael-min.js" ></script> 
      <script type="text/javascript" src="js/jsphylosvg-min.js"></script> 
      
      <script type="text/javascript">
      window.onload = function(){
          var dataObject = { newick: '#{newick_inner_content};' };
          phylocanvas = new Smits.PhyloCanvas(
            dataObject,
            'svgCanvas', 
            #{size_x}, #{size_y},
            'circular'
            }
          );
      };
      </script>
    </head>
    <body>
      <div id="svgCanvas"> </div>
    </body>
    </html>
    RESULT
  end
  
  def create_newick_inner_content(ind = clusterer.root_node, branch_len_meth)
    if clusterer.leaf?(ind)
      compact_name( clusterer.names[ind] )
    else
      ind1,ind2 = clusterer.children(ind)
      dist, dist1, dist2 = clusterer.send(branch_len_meth, ind), clusterer.send(branch_len_meth, ind1), clusterer.send(branch_len_meth, ind2)
      len_1 = (dist - dist1).round(3)
      len_2 = (dist - dist2).round(3)
      "(#{create_newick_inner_content(ind1, branch_len_meth)}:#{len_1},#{create_newick_inner_content(ind2, branch_len_meth)}:#{len_2})"
    end
  end
  
end