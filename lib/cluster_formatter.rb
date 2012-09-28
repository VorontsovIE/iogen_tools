class ClusterFormatter
  attr_accessor :clusterer
  def initialize(clusterer)
    @clusterer = clusterer
  end
  def compact_name(name)
    mtch = /(?<rep>.+).+?\+\k<rep>/.match(name)
    mtch ? name.gsub(/\+#{mtch[:rep]}/,'+') : name
  end
  
end