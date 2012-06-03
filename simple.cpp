
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <inttypes.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <exception>
#include <algorithm>

namespace blossom {

class exception_t: public std::exception {
public:
	exception_t(const char* fmt,...) {
		va_list ap;
		va_start(ap,fmt);
		static char buf[1024];
		*buf = 0;
		vsnprintf(buf,sizeof(buf),fmt,ap);
		va_end(ap);
		_what = std::string(buf);
		//VALGRIND_PRINTF_BACKTRACE("Blossom Exception: %s\n",buf);
	}
	virtual ~exception_t() throw() {}
	virtual const char* what() const throw() { return _what.c_str(); }
protected:
	std::string _what;
};

struct FASTQ_t {
	FASTQ_t(const std::string& n,const std::string& s,const std::string& q): name(n), seq(s), quality(q) {}
	FASTQ_t(const FASTQ_t& copy): name(copy.name), seq(copy.seq), quality(copy.quality) {}
	static FASTQ_t* read(std::istream& in);
	const std::string name, seq, quality;
};

FASTQ_t* FASTQ_t::read(std::istream& in) {
	std::string name, seq, plus, quality;
	std::getline(in,name);
	if(!name.size() || name.at(0) != '@') if(in.eof()) return NULL; else throw exception_t("bad FASTQ name"); 
	std::getline(in,seq);
	if(!seq.size()) throw exception_t("bad FASTQ %s seq length",name.c_str());
	for(int i=0; i<seq.size(); i++)
		switch(seq.at(i)) {
		case 'A': case 'C': case 'G': case 'T': case 'N': break;
		default: throw exception_t("error reading FASTQ base %c at %d in %s",
				seq.at(i),i,name.c_str());
		}
	std::getline(in,plus);
	if(!plus.size() || plus.at(0) != '+') throw exception_t("bad FASTQ %s plus",name.c_str());
	std::getline(in,quality);
	if(quality.size() != seq.size()) throw exception_t("bad FASTQ %s quality length",name.c_str());
	for(int i=0; i<quality.size(); i++)
		if(quality.at(i) < 33 || quality.at(i) > 126) throw exception_t("bad FASTQ %s quality",name.c_str());
	if(!in.good()) throw exception_t("bad FASTQ %s read",name.c_str());
	return new FASTQ_t(name.c_str()+1,seq,quality);
}

class genome_t {
public:
	enum base_t {
		T = 0,
		C = 1,
		A = 2,
		G = 3
	};
	struct seq_t {
		std::string name, data;
	};
	typedef std::vector<seq_t> seqs_t;
	seqs_t seqs;
	typedef std::vector<const char*> sa_t;
	sa_t sa;
	typedef std::vector<uint32_t> cumulative_t;
	cumulative_t cumulative;
	static genome_t* load_FASTA(const std::string& filename);
	virtual ~genome_t() {}
	void match_FASTQs(const std::string& filename) const;
	void match(const FASTQ_t& fq) const;
private:
	genome_t(): prefix_len(10) {}
	int prefix_len;
	void build_tables();
};

genome_t* genome_t::load_FASTA(const std::string& filename) {
	std::auto_ptr<genome_t> ret(new genome_t());
	std::ifstream f(filename.c_str());
	if(!f.is_open()) throw exception_t("could not load FASTA %s",filename.c_str());
	std::clog << "loading FASTA file " << filename << std::endl;
	int line_count;
	std::string line;
	seq_t seq;
	for(line_count=0; f.good(); line_count++) {
		std::getline(f,line);
		if(!line.size()) continue;
		if(!line_count) {
			if(line.at(0) != ';' && line.at(0) != '>')
				throw exception_t("FASTA file did not start with > line");
			seq.name = line.c_str()+1;
		} else if(line.at(0) != ';') {
			for(int i=0; i<line.size(); i++)
				switch(line.at(i)) {
				case 'A': case 'C': case 'G': case 'T': break;
				default: throw exception_t("error reading FASTA base %c at %d:%d in %s",
						line.at(i),line_count,i,filename.c_str());
				}
			seq.data += line;
		}
	}
	if(!f.eof()) throw exception_t("failed to load FASTA %s",filename.c_str());
	ret->seqs.push_back(seq);
	std::clog << "loaded " << line_count << " lines with " << seq.data.size() << " bases from FASTA " << seq.name << std::endl;
	ret->build_tables();
	return ret.release();
}

static bool str_less(const char *c1, const char *c2) { return strcmp(c1, c2) < 0; }

void genome_t::build_tables() {
	sa.clear();
	for(seqs_t::const_iterator s=seqs.begin(); s!=seqs.end(); s++) {
		if(s->data.size() < prefix_len) throw exception_t("sequence %s too short",s->name.c_str());
		sa.reserve(sa.size() + s->data.size() - prefix_len);
		for(const char *i=s->data.c_str(), *const sentinal=s->data.c_str()+s->data.size()-prefix_len; i<sentinal; i++)
			sa.push_back(i);
	}
	std::clog << "building suffix array for " << sa.size() << " entries" << std::endl;
	std::stable_sort(sa.begin(),sa.end(),str_less); // libdivsufsort etc
	std::clog << "building cumulative array for prefix_len=" << prefix_len << std::endl;
	cumulative.clear();
	cumulative.resize(1<<(prefix_len*2),0);
	for(seqs_t::const_iterator s=seqs.begin(); s!=seqs.end(); s++) {
		uint32_t key = 0;
		for(const char *i=s->data.c_str(), *const start=s->data.c_str()+prefix_len, *const sentinal=s->data.c_str()+s->data.size(); i<sentinal; i++) {
			key <<= 2;
			switch(*i) {
			case 'A': key |= A; break;
			case 'C': key |= C; break;
			case 'G': key |= G; break;
			case 'T': key |= T; break;
			default: throw exception_t("bad base %c",*i);
			}
			key &= (1<<(prefix_len*2))-1;
			if(i>start)
				cumulative.at(key)++;
		}
	}
	cumulative_t::value_type sum=0, zero=0, count1=0, count2=0, count4=0, count8=0, big=0, max=0;
	for(cumulative_t::iterator i=cumulative.begin(); i!=cumulative.end(); i++) {
		if(!*i) zero++;
		else if(*i==1) count1++;
		else if(*i==2) count2++;
		else if(*i<=4) count4++;
		else if(*i<=8) count8++;
		else big++;
		max = std::max(max,*i);
		cumulative_t::value_type count = *i;
		*i = sum;
		sum += count;
	}
	cumulative.push_back(sum); // terminator
	std::clog << "cumulative stats: " << cumulative.size() << " zero=" << zero << ", 1=" << count1 <<
		", 2=" << count2 << ", 4=" << count4 << ", 8=" << count8 <<
		", big=" << big << ", max=" << max << ", total=" << cumulative.back() << std::endl;
	std::clog << "RAM: " << (sa.size()*sizeof(sa_t::value_type) + cumulative.size()*sizeof(cumulative_t::value_type)) + sa.size() + prefix_len*seqs.size() << " bytes" << std::endl;
}

void genome_t::match_FASTQs(const std::string& filename) const {
	std::ifstream f(filename.c_str());
	if(!f.is_open()) throw exception_t("could not load FASTQ %s",filename.c_str());
	std::clog << "loading FASTQ file " << filename << std::endl;
	for(std::auto_ptr<FASTQ_t> fq(FASTQ_t::read(f)); fq.get(); fq.reset(FASTQ_t::read(f)))
		match(*fq);
}

void genome_t::match(const FASTQ_t& fq) const {
	std::clog << "matching " << fq.name << std::endl;
	//http://bowtie-bio.sourceforge.net/manual.shtml#the--n-alignment-mode
}

} // namespace blossom

int main(int argc,char** args) {
	using namespace blossom;
	std::auto_ptr<genome_t> e_coli(genome_t::load_FASTA("../bowtie-0.12.8/genomes/NC_008253.fna"));
	e_coli->match_FASTQs("../bowtie-0.12.8/reads/e_coli_1000.fq");
	return 0;
}
