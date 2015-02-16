#
# = bio/appl/sam/report.rb - SAM parser
#
# License:: The Ruby License
#
#  $Id:$
#
#
# == References
# * http://samtools.sourceforge.net/

require 'bio'

module Bio
  class Sam

    # Bio::Sam::Report is a SAM parser class.
    #
    class Report #< DB
      # Delimiter of each entry. Bio::FlatFile uses it.
      # In Bio::Blat::Report, it it nil (1 entry 1 file).
      DELIMITER = RS = nil # 1 file 1 entry

      # Splitter for Bio::FlatFile
      FLATFILE_SPLITTER = Bio::FlatFile::Splitter::LineOriented

      def initialize(text = '')
        flag = false
        @head = []
        @hits = []
        text.each_line do |line|
          if flag then
            @hits << Hit.new(line)
          else
            # for headerless data
            if /^[^@]/ =~ line then
              flag = true
              redo
            end
            @head << line
          end
        end
      end

      # Adds a header line
      # Returns the header_lines
      def add_header_line(line)
        line = line.chomp
        if line =~ /^@/
          @header_lines ||= []
          @header_lines.push line
          return self
        else
          return false
        end
      end

      # Adds a line to the entry if the given line 
      # is derived from the same template of the current entry.
      # If the current entry (self) is empty, or the line has the same
      # template name, the line is added and returns self.
      # Otherwise, returns false (the line is not added).
      def add_line(line)
        if /\A\s*\z/ =~ line then
          return @hits.empty? ? self : false
        end
        hit = Hit.new(line.chomp)
        if @hits.empty? or @hits.first.template_id == hit.template_id then
          @hits.push hit
          return self
        else
          return false
        end
      end

      # hits of the result.
      # Returns an Array of Bio::Sam::Report::Hit objects.
      attr_reader :hits

      # Hit class for the Sam parser.
      class Hit
        # Creates a new Hit object from a SAM record.
        # It is designed to be called internally from Bio::Sam::Report object.
        # Users shall not use it directly.
        def initialize(str)
          @data = str.chomp.split(/\t/)
        end

        # Raw data of the hit.
        attr_reader :data

        # Returns strand information of the hit.
        # Returns '+' or '-'.
        # This would be a Bio::Blat specific method.
        def strand;      (@data[1].to_i & 0x10 != 0) ? '-' : '+';       end
        def template_id; @data[0]; end
        alias query_id template_id
        def target_id; @data[2]; end
        def flag; @data[1].to_i; end
        def pos; @data[3].to_i;end
        def mapq; @data[4].to_i; end
        def cigar; @data[5];end
        def rnext; @data[6]; end
        def pnext; @data[7].to_i; end
        def tlen; @data[8].to_i; end
        def seq; @data[9]; end
        def qual_str; @data[10]; end

      end #class Hit
    end #class Report
  end #class Sam
end #module Bio

=begin

= Bio::Sam::Report

  SAM parser.

= References

 * http://samtools.sourceforge.net/
=end
