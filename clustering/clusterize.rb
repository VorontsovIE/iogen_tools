# ruby clusterize.rb [--log log_file.log] <matrix-txt-file> <motif-names-yaml-file> <output-folder> [cluster-yaml-dump]
# cluster yaml dump used in a pair of ways:
#  - when specified dump exists - it's loaded (in such case time-consuming clusterization stage's eliminated)
#  - when dump not exists - built clusterer'll be dumped to specified file (so that next time clusterization run, it could immediately get clusterization tree)

require_relative 'lib/clusterer'
require 'yaml'
require 'fileutils'
require 'optparse'


# Returns sorted list of cutoffs with number of clusters at this cutoff and number of annotated clusters
# Format: [[cutoff_1, num_clusters_1, num_annotated_clusters_1], ...]
# Selection condition should be specified with a block (because there are different known-denovo recognition methods)
# Block gets cluster (names of motifs in cluster) and should return true iff cluster is annotated
def calculate_statistics_by_possible_cutoffs(clusterer, criterium, &block)
  dists = (0..clusterer.root_node).map{|ind| clusterer.send(criterium,ind)}.uniq.sort
  result = []
  dists.each{|cutoff|
    clusters = clusterer.get_clusters_names(&clusterer.cutoff_criterium(criterium, cutoff))
    annotated_clusters = clusters.select(&block)
    result << [cutoff, clusters.size, annotated_clusters.size]
  }
  File.open("#{criterium}_cutoffs.txt",'w'){|f| f.puts result.map{|cutoff_stat| cutoff_stat.join("\t")}}
  result
end

# Looks for a cutoff that yields maximal number of annotated clusters (consisting of both known and denovo motif).
# Returns least possible cutoff (prefer more compact clusters)
def best_cutoff(clusterer, criterium, &block)
  statistics = calculate_statistics_by_possible_cutoffs(clusterer, criterium, &block)
  max_num_annotated_clusters = statistics.map{|cutoff, num_clusters, num_annotated_clusters| num_annotated_clusters }.max
  
  statistics.find{|cutoff, num_clusters, num_annotated_clusters| num_annotated_clusters == max_num_annotated_clusters }.tap do |x|
    puts "best #{criterium} cutoff: #{x[0]} -- num clusters: #{x[1]} -- num annotated_clusters: #{x[2]}"
  end.first
end

options = { }
OptionParser.new{|cmd|
  cmd.on('-l', '--log LOG_FILE', 'log-file of clusterization process (by default stderr used)'){ |log_file|
    options[:log_file] = log_file
  }
}.parse!

matrix_filename = ARGV.shift        # distance_matrix/distance_macroape.txt
names_filename = ARGV.shift         # distance_matrix/motifs_order.yaml
output_folder = ARGV.shift          # distance_matrix/clustering_results
cluster_dump_filename = ARGV.shift  # distance_matrix/cluster.yaml

raise 'matrix filename not specified'  unless matrix_filename
raise "matrix file #{matrix_filename} not exist"  unless File.exist?(matrix_filename)
raise 'names filename not specified'  unless names_filename
raise "names file #{names_filename} not exist"  unless File.exist?(names_filename)
raise 'output folder not specified'  unless output_folder

FileUtils.mkdir_p(output_folder)  unless Dir.exist? output_folder
FileUtils.mkdir_p(File.dirname(cluster_dump_filename))  if cluster_dump_filename && ! Dir.exist?(File.dirname(cluster_dump_filename))

names = YAML.load_file(names_filename)
distance_matrix = load_matrix_from_file(matrix_filename)


if cluster_dump_filename
  if File.exist?(cluster_dump_filename)
    clusterer = Clusterer.load(distance_matrix, cluster_dump_filename)
  else
    clusterer = Clusterer.new(distance_matrix, :average_linkage, names)
    clusterer.logger = Logger.new(options[:log_file])  if options[:log_file]
    clusterer.make_linkage
    clusterer.dump(cluster_dump_filename)
  end
else
  clusterer = Clusterer.new(distance_matrix, :average_linkage, names)
  clusterer.make_linkage
end

newick_formatter = ClusterNewickFormatter.new(clusterer, :link_length)
xml_formatter = ClusterXMLFormatter.new(clusterer, :link_length, 0.1)
File.open("#{output_folder}/macroape_linklength.html",'w'){|f| f << newick_formatter.create_newick_html()}
File.open("#{output_folder}/macroape_linklength.newick",'w'){|f| f << newick_formatter.content()}

File.open("#{output_folder}/macroape_linklength.xml",'w'){|f| f << xml_formatter.content()}
File.open("#{output_folder}/macroape_linklength_w_xml.html",'w'){|f| f << xml_formatter.create_html_connected_to_xml("#{output_folder}/macroape_linklength.xml") }


distance_macroape_cutoff_grid = (0.90...1).step(0.01).to_a

distance_macroape_cutoff_grid << best_cutoff(clusterer, :subtree_max_distance){|cluster| cluster.any?{|name| name =~ /KNOWN/ }  &&  cluster.any?{|name| name !~ /KNOWN/ }  }
distance_macroape_cutoff_grid << best_cutoff(clusterer, :link_length){|cluster| cluster.any?{|name| name =~ /KNOWN/ }  &&  cluster.any?{|name| name !~ /KNOWN/ }  }


clusters_macroape_linklength, clusters_macroape_maxdist = {}, {}
distance_macroape_cutoff_grid.each { |cutoff|
  clusters_macroape_linklength[cutoff] = clusterer.get_clusters_names(&clusterer.cutoff_criterium(:link_length, cutoff))
  clusters_macroape_maxdist[cutoff] = clusterer.get_clusters_names(&clusterer.cutoff_criterium(:subtree_max_distance, cutoff))
}

clusters_macroape_linklength.each do |cutoff, clusters|
  File.open("#{output_folder}/macroape_linklength_cluster_names (#{cutoff.round(4)} - #{clusters.size}).txt",'w'){|f|
    clusters.each{|clust| f.puts clust.join("\t") }
  }    
end

clusters_macroape_maxdist.each do |cutoff, clusters|
  File.open("#{output_folder}/macroape_maxdist_cluster_names (#{cutoff.round(4)} - #{clusters.size}).txt",'w'){|f| 
    clusters.each{|clust| f.puts clust.join("\t")}
  }
end
