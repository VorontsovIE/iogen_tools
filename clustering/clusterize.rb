# ruby clusterize.rb [--log log_file.log] <matrix-txt-file> <motif-names-yaml-file> <output-folder> [cluster-yaml-dump]
# cluster yaml dump used in a pair of ways:
#  - when specified dump exists - it's loaded (in such case time-consuming clusterization stage's eliminated)
#  - when dump not exists - built clusterer'll be dumped to specified file (so that next time clusterization run, it could immediately get clusterization tree)

require_relative 'lib/clusterer'
require 'yaml'
require 'fileutils'
require 'optparse'

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

options = { }
OptionParser.new{|cmd|
  cmd.on('-l', '--log LOG_FILE', 'log-file of clusterization process (by default stderr used)'){ |log_file|
    options[:log_file] = log_file
  }
}.parse!

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
end

# File.open("#{output_folder}/macroape_linklength.newick",'w'){|f| f << clusterer.create_newick_inner_content(:link_length)}

distance_macroape_cutoff_grid = (0.90...1).step(0.01).to_a

clusters_macroape_linklength, clusters_macroape_maxdist = {}, {}
distance_macroape_cutoff_grid.each { |cutoff|
  clusters_macroape_linklength[cutoff] = clusterer.get_clusters_names(&clusterer.cutoff_criterium(:link_length, cutoff))
  clusters_macroape_maxdist[cutoff] = clusterer.get_clusters_names(&clusterer.cutoff_criterium(:subtree_max_distance, cutoff))
}

clusters_macroape_linklength.each do |cutoff, clusters|
  File.open("#{output_folder}/macroape_linklength_cluster_names (#{cutoff.round(3)} - #{clusters.size}).txt",'w'){|f|
    clusters.each{|clust| f.puts clust.join("\t") }
  }    
end

clusters_macroape_maxdist.each do |cutoff, clusters|
  File.open("#{output_folder}/macroape_maxdist_cluster_names (#{cutoff.round(3)} - #{clusters.size}).txt",'w'){|f| 
    clusters.each{|clust| f.puts clust.join("\t")}
  }
end
