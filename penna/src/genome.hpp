#ifndef GENOME_HPP
#define GENOME_HPP

#include <bitset>

namespace Penna{

using age_t = unsigned int;

class Genome{
	
	public:
		static age_t const genome_size = 64;
		using cont_gen = std::bitset<genome_size>;

		Genome() {
			for(std::size_t k = 0; k < genome_size; ++k){
				genes_[k] = false;
			}
		}

		static void  set_mutation_rate(age_t);
		static age_t get_mutation_rate();

		cont_gen get_genome() const;

		void mutate();

		age_t count_bad(age_t) const;

	private:
		static age_t mutation_rate_;
		cont_gen genes_;
};

} // Penna namespace

#endif // _GENOME_HPP_