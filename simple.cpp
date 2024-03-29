
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
		A, C, G, T, // alphabetic so strcmp works too
	};
	struct seq_t {
		std::string name;
		size_t ofs, len;
	};
	typedef std::vector<seq_t> seqs_t;
	seqs_t seqs;
	std::string data;
	typedef std::vector<uint32_t> sa_t;
	sa_t sa;
	typedef std::vector<uint32_t> cumulative_t;
	cumulative_t cumulative;
	struct match_t {
		match_t(): ofs(0), len(0), score(0) {}
		match_t(sa_t::value_type o,uint32_t l,uint32_t s): ofs(o), len(l), score(s) {}
		sa_t::value_type ofs;
		uint32_t len;
		uint32_t score;
	};
	typedef std::vector<match_t> matches_t;
	static genome_t* load_FASTA(const std::string& filename);
	static genome_t* load_blossom(const std::string& filename);
	void save_to_blossom(const std::string& filename) const;
	void load_from_blossom(const std::string& filename);
	virtual ~genome_t() {}
	void match_FASTQs(const std::string& filename) const;
	void match(const FASTQ_t& fq,matches_t& matches) const;
private:
	genome_t(): prefix_len(10) {}
	int prefix_len;
	void build_tables();
};

genome_t* genome_t::load_blossom(const std::string& filename) {
	std::auto_ptr<genome_t> ret(new genome_t());
	ret->load_from_blossom(filename);
	return ret.release();
}

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
			seq.ofs = ret->data.size();
			seq.len = 0;
		} else if(line.at(0) != ';') {
			for(int i=0; i<line.size(); i++)
				switch(line.at(i)) {
				case 'A': case 'C': case 'G': case 'T': break;
				default: throw exception_t("error reading FASTA base %c at %d:%d in %s",
						line.at(i),line_count,i,filename.c_str());
				}
			seq.len += line.size();
			ret->data += line;
		}
	}
	if(!f.eof()) throw exception_t("failed to load FASTA %s",filename.c_str());
	ret->seqs.push_back(seq);
	std::clog << "loaded " << line_count << " lines with " << seq.len << " bases from FASTA " << seq.name << std::endl;
	ret->build_tables();
	return ret.release();
}

template<typename T> void write(std::ostream& f,const T& t) { f.write((const char*)&t,sizeof(T)); }
template<typename T> T read(std::istream& f) { T t; f.read((char*)&t,sizeof(T)); return t; }
static const char blossom_signature_text[] = {'b','l','o','s','s','o','m',0};
static const uint32_t blossom_signature_endian = 0x01F0;

void genome_t::save_to_blossom(const std::string& filename) const {
	std::ofstream f(filename.c_str());
	if(!f.is_open()) throw exception_t("could not create blossom file %s",filename.c_str());
	std::clog << "saving blossom file " << filename << std::endl;
	write(f,blossom_signature_text);
	write(f,blossom_signature_endian);
	write<uint8_t>(f,1); // version
	write<uint8_t>(f,sizeof(sa_t::value_type));
	write<uint8_t>(f,sizeof(cumulative_t::value_type));
	write<uint32_t>(f,seqs.size());
	for(seqs_t::const_iterator s=seqs.begin(); s!=seqs.end(); s++) {
		write<uint32_t>(f,s->name.size());
		f.write(s->name.c_str(),s->name.size());
		write<uint32_t>(f,s->ofs);
		write<uint32_t>(f,s->len);
	}
	write<uint32_t>(f,data.size());
	f.write(data.c_str(),data.size());
	write<char>(f,prefix_len);
	write<uint32_t>(f,sa.size());
	f.write((const char*)&sa.front(),sa.size()*sizeof(sa_t::value_type));
	write<uint32_t>(f,cumulative.size());
	f.write((const char*)&cumulative.front(),cumulative.size()*sizeof(cumulative_t::value_type));
}

void genome_t::load_from_blossom(const std::string& filename) {
	std::ifstream f(filename.c_str());
	if(!f.is_open()) throw exception_t("could not open blossom file %s",filename.c_str());
	std::clog << "loading blossom file " << filename << std::endl;
	char sig_check[sizeof(blossom_signature_text)];
	f.read(sig_check,sizeof(blossom_signature_text));
	if(memcmp(sig_check,blossom_signature_text,sizeof(blossom_signature_text))) throw exception_t("not a blossom file");
	if(read<uint32_t>(f) != blossom_signature_endian) throw exception_t("bad blossom endian");
	if(read<uint8_t>(f) != 1) throw exception_t("bad blossom file version");
	if(read<uint8_t>(f) != sizeof(sa_t::value_type)) throw exception_t("incompatible blossom SA type");
	if(read<uint8_t>(f) != sizeof(cumulative_t::value_type)) throw exception_t("incompatible blossom cumulative counter type");
	// at this point everything seems set up, so go blank everything ready
	seqs.clear();
	sa.clear();
	cumulative.clear();
	data.clear();
	// and read it in
	for(uint32_t num_seqs = read<uint32_t>(f); num_seqs--;) {
		seq_t seq;
		seq.name.resize(read<uint32_t>(f));
		f.read(&seq.name.at(0),seq.name.size());
		seq.ofs = read<uint32_t>(f);
		seq.len = read<uint32_t>(f);
		seqs.push_back(seq);
	}
	data.resize(read<uint32_t>(f));
	f.read(&data.at(0),data.size());
	prefix_len = read<char>(f);
	sa.resize(read<uint32_t>(f));
	f.read((char*)&sa.at(0),sa.size()*sizeof(sa_t::value_type));
	cumulative.resize(read<uint32_t>(f));
	f.read((char*)&cumulative.at(0),cumulative.size()*sizeof(cumulative_t::value_type));
	std::clog << data.size() << " bytes of data, " << sa.size() << " SA entries and " << cumulative.size() << " entries in the cumulative lookup" << std::endl;
}

struct str_less {
	str_less(const char* p_): p(p_) {}
	bool operator()(const genome_t::sa_t::value_type a, const genome_t::sa_t::value_type b) { return strcmp(p+a,p+b) < 0; }
	const char * const p;
};

void genome_t::build_tables() {
	sa.clear();
	for(seqs_t::const_iterator s=seqs.begin(); s!=seqs.end(); s++) {
		if(s->len < prefix_len) throw exception_t("sequence %s too short",s->name.c_str());
		sa.reserve(sa.size() + s->len - prefix_len);
		for(sa_t::value_type i=s->ofs, sentinal=s->ofs+s->len-prefix_len; i<sentinal; i++)
			sa.push_back(i);
	}
	std::clog << "building suffix array for " << sa.size() << " entries" << std::endl;
	std::stable_sort(sa.begin(),sa.end(),str_less(data.c_str())); // libdivsufsort etc
	for(sa_t::iterator i=sa.begin(); i!=sa.end(); i++) // add offset; does this make it a prefix array? :)
		*i += prefix_len;
	std::clog << "building cumulative array for prefix_len=" << prefix_len << std::endl;
	cumulative.clear();
	cumulative.resize(1<<(prefix_len*2),0);
	for(seqs_t::const_iterator s=seqs.begin(); s!=seqs.end(); s++) {
		uint32_t key = 0;
		for(const char *i=data.c_str()+s->ofs, *const start=data.c_str()+s->ofs+prefix_len, *const sentinal=data.c_str()+s->ofs+s->len; i<sentinal; i++) {
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
	int exacts = 0;
	matches_t matches;
	for(std::auto_ptr<FASTQ_t> fq(FASTQ_t::read(f)); fq.get(); fq.reset(FASTQ_t::read(f))) {
		matches.clear();
		match(*fq,matches);
		bool exact = false;
		for(matches_t::const_iterator i=matches.begin(); i!=matches.end(); i++)
			if(i->len==fq->seq.size()) {
				std::clog << fq->name << '\t' << i->ofs << std::endl;
				exact = true;
			}
		if(exact) exacts++;
	}
	std::clog << exacts << " exact matches" << std::endl;
}

void genome_t::match(const FASTQ_t& fq,matches_t& matches) const {
	//http://bowtie-bio.sourceforge.net/manual.shtml#the--n-alignment-mode
	// exact match
	uint32_t key = 0;
	for(int i=0; i<prefix_len; i++) {
		key <<= 2;
		switch(fq.seq.at(i)) {
		case 'A': key |= A; break;
		case 'C': key |= C; break;
		case 'G': key |= G; break;
		case 'T': key |= T; break;
		case 'N': {
			//std::clog << "abandoned unknown" << std::endl;
			return;
		} break;
		default: throw exception_t("bad base %c at FASTQ %d",fq.seq.at(i),i);
		}
		key &= (1<<(prefix_len*2))-1;
	}
	const cumulative_t::value_type sa_start = cumulative.at(key), sa_stop = cumulative.at(key+1);
	int best = 0;
	cumulative_t::value_type best_ofs;
	for(cumulative_t::value_type sa_ofs = sa_start; sa_ofs < sa_stop; sa_ofs++) {
		__builtin_prefetch(data.c_str()+sa.at(sa_ofs+1)+prefix_len);
		int match = prefix_len;
		for(const char* a=fq.seq.c_str()+prefix_len, *b=data.c_str()+sa.at(sa_ofs); *a==*b && *a; a++,b++,match++);
		if(match > best) {
			best = match;
			best_ofs = sa.at(sa_ofs)-prefix_len;
		}
	}
	if(best)
		matches.push_back(match_t(best_ofs,best,best));
}

} // namespace blossom

int main(int argc,char** args) {
	std::ios_base::sync_with_stdio(false);	
	using namespace blossom;
	std::auto_ptr<genome_t> e_coli;
	if(false) {
		e_coli.reset(genome_t::load_FASTA("../bowtie-0.12.8/genomes/NC_008253.fna"));
		e_coli->save_to_blossom("e_coli.blossom");
	} else
		e_coli.reset(genome_t::load_blossom("e_coli.blossom"));
	e_coli->match_FASTQs("../bowtie-0.12.8/reads/e_coli_1000.fq");
	return 0;
}
