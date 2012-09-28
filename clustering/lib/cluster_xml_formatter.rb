require_relative 'cluster_formatter'

class ClusterXMLFormatter < ClusterFormatter
  def create_html_connected_to_xml(xml_file, size={})
    size_x = size[:x] || 2000
    size_y = size[:y] || 2000
  
    <<-RESULT
    <html>
    <head>
      <link type="text/css" rel="stylesheet" href="js/yui/build/cssfonts/fonts-min.css" /> 
      <script type="text/javascript" src="js/raphael/raphael-min.js" ></script> 
      <script type="text/javascript" src="js/yui/build/yui/yui.js"></script> 
      <script type="text/javascript" src="js/jsphylosvg-min.js"></script> 
      <!-- unitip -->
      <link rel="stylesheet" type="text/css" href="js/unitip/css/unitip.css" > 
      <script type="text/javascript" src="js/unitip/js/unitip.js"></script> 	
      <script type="text/javascript">
      window.onload = function(){
        YUI().use('oop', 'json-stringify', 'io-base', 'event', 'event-delegate', function(Y){
          var uri = "#{xml_file}";
          function complete(id, o, args) {
            var data = o.responseXML; // Response data.
            var dataObject = {
                  xml: data,
                  fileSource: true
                };
            phylocanvas = new Smits.PhyloCanvas(
              dataObject,
              'svgCanvas', 
              #{size_x}, #{size_y},
              'circular'
            );
            
            document.getElementById('svgCode').value = phylocanvas.getSvgSource();
            init(); //unitip
          };
          Y.on('io:complete', complete, Y);
          var request = Y.io(uri);
        });
      };
      </script>
    </head>
    <body>
      <div id="svgCanvas"> </div>
      <textarea id="svgCode"></textarea>
    </body>
    </html>
    RESULT
  end
  
  def create_xml(branch_len_meth, cutoff)
    clusters = clusterer.subtree_clusters(&clusterer.cutoff_criterium(branch_len_meth, cutoff))
    xml_inner_content = create_xml_inner_content(branch_len_meth, cutoff, clusters)
    
    <<-RESULT
    <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
    <phylogeny rooted="false">
    <render>
      <!-- <charts>
        <component type="binary" thickness="10" isInternal="true" />
      </charts> -->
      <styles>
        <first_cluster_group fill='#9F9' stroke='#FFF' />
        <second_cluster_group fill='#F99' stroke='#FFF' />
      </styles>
      <parameters>
        <circular>
            <innerCircleRadius>100</innerCircleRadius>
        </circular>
      </parameters>
    </render>
    <clade>
    <branch_length> 0.0 </branch_length>
    <clade>
    <branch_length> 0.0 </branch_length>
    #{ xml_inner_content }
    </clade>
    </clade>
    </phylogeny>
    </phyloxml>
    RESULT
  end
  
  # block: ind in linkage tree  -->  dist 
  def create_xml_inner_content(ind = clusterer.root_node, branch_len_meth, cutoff, clusters)
    if clusterer.leaf?(ind)
      stylename = clusters.find_index{|cluster| cluster.include? ind} % 2 == 0 ? 'first_cluster_group' : 'second_cluster_group'
      name = clusterer.names[ind]
      <<-RETSTRING
        <name bgStyle="#{stylename}">#{compact_name(name)}</name>
        <annotation>
        <desc>&lt;img src=\"http://autosome.ru/hocomoco/hocomoco_f/motifs/#{name}_thumb.jpg\"/&gt;</desc>
        <uri>http://autosome.ru/hocomoco/hocomoco_ad/motif_details/#{name}.html</uri>
        </annotation>
        <chart>
          <component>#{stylename}</component>
        </chart>
      RETSTRING
    else
      ind1, ind2 = clusterer.children(ind)
      dist, dist1, dist2 = clusterer.send(branch_len_meth, ind), clusterer.send(branch_len_meth, ind1), clusterer.send(branch_len_meth, ind2)
      len_1 = (dist - dist1).round(3)
      len_2 = (dist - dist2).round(3)
      <<-RETSTRING
        <clade>
        <branch_length>#{len_1}</branch_length>
        #{create_xml_inner_content(ind1, branch_len_meth, cutoff, clusters)}
        </clade>
        <clade>
        <branch_length>#{len_2}</branch_length>
        #{create_xml_inner_content(ind2, branch_len_meth, cutoff, clusters)}
        </clade>
      RETSTRING
    end
  end
end