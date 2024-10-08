#!/usr/bin/env ruby

require 'bio'
require 'optparse'
require_relative 'samreport'

$start_pos_limit=0
$report_threshold= 3
$max_n = 3
$filter_t1l = 1
$filter_t1r = 0
$filter_t2l = 1
$filter_t2r = 0
OptionParser.new do |opt|
  opt.on('-r FASTAFILE', '--ref=FILENAME'){|fastafile| 
    $reffastafile=fastafile
  }#'PhiX.fa'
  opt.on('-s SAMFILE', '--sam=FILENAME'){|samfile| 
    $samfile=samfile
  }#'PhiX1Mpebwa.sam'
  opt.on('-q MIN_QV', '--min-qv=INTEGER'){|v| 
    $threshold = v.to_i
  }
  opt.on('-t REPORT_THRESHOLD', '--report-threshold=INTEGER'){|v| 
    $report_threshold = v.to_i
  }
  opt.on('-n MAX_N', '--max-N=INTEGER'){|v| 
    $max_n = v.to_i
  }
  opt.on('-m MAX_FRAG', '--max-fragment-size=INTEGER'){|lb|
    $length_barrier=lb.to_i
  }
  opt.on('-p START_POS', '--start-pos=INTEGER'){|sp|
    $start_pos_limit=sp.to_i
# if the match begins before this point the hit should be discarded
# start position of the left primer should be set
  }
  opt.on('-e END_POS', '--end-pos=INTEGER'){|ep|
    $end_pos_limit=ep.to_i
  }
# if the match extends beyond this point the hit should be discarded
  opt.on('-c START_POS_CUT', '--left-primer-end=INTEGER'){|sc|
    $start_pos_cut=sc.to_i
  }
# if the match begins prior to this point the region up to this point
# should be ignored as this is priming site
# the end position of left primer of the amplicon
# 0 if no pcr took place
  opt.on('-d END_POS_CUT', '--right-primer-end=INTEGER'){|ep|
    $end_pos_cut = ep.to_i
  }
  opt.on('-g 3D-PCR end 1', '--exclude-target1-l=INTEGER'){|f|
    $filter_t1l = f.to_i
  } 
  opt.on('-h 3D-PCR end 2', '--exclude-target1-r=INTEGER'){|f|
    $filter_t1r = f.to_i
  } 
  opt.on('-i 3D-PCR end 1', '--exclude-target2-l=INTEGER'){|f|
    $filter_t2l = f.to_i
  } 
  opt.on('-j 3D-PCR end 2', '--exclude-target2-r=INTEGER'){|f|
    $filter_t2r = f.to_i
  } 
  opt.parse!(ARGV)
end
puts "Will report reads with more than #{$report_threshold} CT or GA mutations"

def complement(na)
  na.tr("acgt", "tgca")
end

class DiffAccumulator
  def initialize(threshold)
    @threshold = threshold
    @acc = Hash.new
    @acc['a'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    @acc['c'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    @acc['g'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    @acc['t'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    @acc['n'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    @refcontext_ga = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0, nil => 0}
    @refcontext_g = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0, nil => 0}
    @refcontext_ct = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0, nil => 0}
    @refcontext_c = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0, nil => 0}
    @nsigmm = Array.new
    @bsigmm = Array.new
#    @nsigmm['ga'] = Array.new
#    @nsigmm['ct'] = Array.new
#    @nsigmm['gact'] = Array.new
  end
  def updatensigmm(pat, count, ngoodbases)
    return if ngoodbases == 0
    if @nsigmm[count] != nil
      @nsigmm[count] += 1
      @bsigmm[count] += ngoodbases
    else
      @nsigmm[count] = 1
      @bsigmm[count] = ngoodbases
    end
  end
  def count(qseq, qqual, refseq)
    qseq = qseq.downcase
#    puts qseq 
#    puts refseq
#    puts qqual.join(' ')
#    exit
    dif= Hash.new
    dif['a'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    dif['c'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    dif['g'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    dif['t'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    dif['n'] = {'a' => 0, 'c' => 0, 'g' => 0, 't' => 0, 'n' =>0}
    ngoodbases = 0
    (0...[qseq.length, refseq.length].min).each do |i|
      if qqual[i] > @threshold
        ngoodbases += 1
        dif[refseq[i]][qseq[i]] += 1 
        @acc[refseq[i]][qseq[i]] += 1 
        if refseq[i] == 'g' 
          @refcontext_g[refseq[i+1]] += 1 unless i + 1 == qseq.length
          if qseq[i] == 'a'
            @refcontext_ga[refseq[i+1]] += 1 unless i + 1 == qseq.length
          end
        elsif refseq[i] == 'c' 
          @refcontext_c[refseq[i-1]] += 1 unless i == 0
          if qseq[i] == 't'
            @refcontext_ct[refseq[i-1]] += 1 unless i == 0
          end
        end
      end
    end
    max = [dif['g']['a'], dif['c']['t']].max
    updatensigmm('gact', max, ngoodbases)
#    updatensigmm('ga', dif['g']['a'])
#    updatensigmm('ct', dif['c']['t'])
    return [max, dif['g']['a'], dif['c']['t']]
  end
  def get_summary
    @acc
  end
  def reportnsigmm(pat)
    puts "number of read fragment with high quality (QV>#{@threshold}) #{pat} substitution\n"
    puts "#subst\t#reads\t#hqbases"
    ary = @nsigmm
    ary.each_index do |i|
      if ary[i] != nil && ary[i] > 0
        puts "#{i}\t#{ary[i]}\t#{@bsigmm[i]}"
      end
    end
    puts "the nucleotide before the C->T substitution"
    puts "for G->A substitution the complement of the next base is accumulated"
    puts "reference context_ga #{@refcontext_ga}"
    puts "reference context_g #{@refcontext_g}"
    puts "reference context_ct #{@refcontext_ct}"
    puts "reference context_c #{@refcontext_c}"
    puts "both strand combined"
    @refcontext_ct_c = @refcontext_ct.dup
    @refcontext_ga.each_pair{|key,val| 
      next if key == nil
      @refcontext_ct_c[complement(key)] += val
    }
    puts "reference context_ct_c #{@refcontext_ct_c}"
    @refcontext_c_c = @refcontext_c.dup
    @refcontext_g.each_pair{|key,val| 
      next if key == nil
      @refcontext_c_c[complement(key)] += val
    }
    puts "reference context_c_c #{@refcontext_c_c}"
    
  end
end

def makeconsensus(seq1, qual1, seq2, qual2, pdiff, tlen)
  return makeconsensus(seq2, qual2, seq1, qual1, -pdiff, tlen) if pdiff < 0
  seq = "N"*tlen
  qual = Array.new(tlen,2)
  seq[pdiff ... pdiff + seq1.length ] = seq1
  qual[pdiff ... pdiff + qual1.length] = qual1
  (0...seq2.length).each do |i|
    if (seq[i]  == seq2[i]) then
      qual[i] += qual2[i] if qual[i] >2
    elsif (seq[i] == "N") then
      seq[i] = seq2[i]
      qual[i] = qual2[i]
    elsif (seq2[i] == "N") then
      # don't change
    else 
      seq[i] = "N"
      qual[i] = 2
    end
  end
  return SeqQual.new(seq,qual)
end
def clip(seq, cigar)
  cigar =~ /(([0-9]+)S)?([0-9]+)M([0-9]+S)?/
  preskip = $2.to_i
  length = $3.to_i
  seq[0+preskip, length]
end
def clipq(seq, cigar, pos, start_cut, end_cut)
  cigar =~ /(([0-9]+)S)?([0-9]+)M([0-9]+S)?/
  preskip = $2.to_i
  length = $3.to_i
  qs = seq[0+preskip, length].dup
  mlen = start_cut - pos + 1
  if mlen > 0
    qs[0...mlen] = '!' * mlen
  end
  mlen = pos + length - end_cut 
  if mlen > 0
    if mlen < length
      qs[-mlen .. -1] = '!' * mlen
    else
      qs = '!' * length
    end
  end
  qs
end
class Refhandler
  def initialize(path)
    @seqs=Hash.new
    @max_size=0
    Bio::FlatFile.open(nil, path).each do |fe|
      @seqs[fe.entry_id]=fe.naseq
      @max_size = fe.naseq.size if fe.naseq.size > @max_size
    end
  end
  def get_ref(name, start, length)
    @seqs[name].subseq(start, start + length -1)
  end
  attr_reader :max_size
end
class SeqQual
  def initialize(seq,qual)
    @seq = seq
    @qual = qual
  end
  attr_reader :seq, :qual
end
refs = Refhandler.new($reffastafile)
$end_pos_limit = refs.max_size if $end_pos_limit == nil
dc = DiffAccumulator.new($threshold)
i = 0
ff = Bio::FlatFile.open(Bio::Sam::Report, $samfile)
# Main loop
ff.each do |fe|
  i+=1
#  puts i.to_s
#  puts fe.hits.size
  next if fe.hits.size != 2
#  puts fe.hits
  next if fe.hits[0].target_id == "*" or fe.hits[1].target_id == "*"
  next if fe.hits[0].tlen > $length_barrier or 
    fe.hits[0].tlen < -$length_barrier
  next if fe.hits[0].tlen == 0 
  next if fe.hits[0].tlen + fe.hits[1].tlen != 0
  next if fe.hits[0].pos < $start_pos_limit
  next if fe.hits[1].pos < $start_pos_limit
#  puts fe.hits[0].tlen
#  puts fe.hits[1].tlen
  tlen0 = fe.hits[0].tlen
  tlen0 = -tlen0 if tlen0 <0
  tlen1 = fe.hits[1].tlen
  tlen1 = -tlen1 if tlen1 <0
  tlen = tlen1 > tlen0 ? tlen1 : tlen0
  next if fe.hits[0].pos >= $filter_t1l && 
          fe.hits[0].pos + tlen0 - 1 <= $filter_t1r
  next if fe.hits[1].pos >= $filter_t1l &&
          fe.hits[1].pos + tlen1 - 1 <= $filter_t1r
  next if fe.hits[0].pos >= $filter_t2l && 
          fe.hits[0].pos + tlen0 - 1 <= $filter_t2r
  next if fe.hits[1].pos >= $filter_t2l &&
          fe.hits[1].pos + tlen1 - 1 <= $filter_t2r
  next if fe.hits[0].pos + tlen > $end_pos_limit
  pdiff = fe.hits[0].pos - fe.hits[1].pos
  clippedseq1 = clip(fe.hits[0].seq, fe.hits[0].cigar)
  clippedseq2 = clip(fe.hits[1].seq, fe.hits[1].cigar)
  if pdiff >  0
    tlen = pdiff + clippedseq1.length if tlen < pdiff + clippedseq1.length
    tlen = clippedseq2.length if tlen < clippedseq2.length
  else
    tlen = clippedseq1.length if tlen < clippedseq1.length
    tlen = -pdiff + clippedseq2.length if tlen < -pdiff + clippedseq2.length
  end
  pmin = pdiff > 0 ? fe.hits[1].pos : fe.hits[0].pos
  abspdiff = pdiff > 0 ? pdiff : -pdiff
  refseq = refs.get_ref(fe.hits[0].target_id, pmin, tlen)
#  puts pdiff
  clippedqual1 = clipq(fe.hits[0].qual_str, fe.hits[0].cigar, fe.hits[0].pos, $start_pos_cut, $end_pos_cut)
  clippedqual2 = clipq(fe.hits[1].qual_str, fe.hits[1].cigar, fe.hits[1].pos, $start_pos_cut, $end_pos_cut)
  cseq=makeconsensus(clippedseq1, clippedqual1.unpack('c*').map{|x| x - 33}, 
    clippedseq2, clippedqual2.unpack('c*').map{|x| x - 33}, pdiff, tlen)
  next if cseq.seq.count("N") > $max_n
#  p cseq.seq.size, cseq.qual.size
  mmc = dc.count(cseq.seq, cseq.qual, refseq) 
  if (mmc[0] > $report_threshold)
    puts mmc
    puts fe.hits[0].template_id
    puts "#{fe.hits[0].target_id}:#{pmin}-"
    puts refseq
    puts fe.hits[0].pos
    puts fe.hits[0].cigar + fe.hits[0].strand
    print " " * abspdiff if pdiff > 0
    puts clippedseq1
    puts clippedqual1.unpack('c*').map{|x| x - 33}.join(' ')
    puts fe.hits[1].pos
    puts fe.hits[1].cigar + fe.hits[1].strand
    print " " * abspdiff if pdiff  <0
    puts clippedseq2
    puts clippedqual2.unpack('c*').map{|x| x - 33}.join(' ')
    puts "read consensus"
    puts cseq.seq
    puts cseq.qual.join(' ')
    if(mmc[1] > $report_threshold)
      puts "GAdiff #{mmc[1]} #{fe.hits[0].template_id}"
    end
    if(mmc[2] > $report_threshold)
      puts "CTdiff #{mmc[2]} #{fe.hits[0].template_id}"
    end
  end
end
puts dc.get_summary
dc.reportnsigmm('gact')
#dc.reportnsigmm('ga')
#dc.reportnsigmm('ct')
