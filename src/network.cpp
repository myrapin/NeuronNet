#include "network.h"
#include "random.h"

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}



std::vector<std::pair<size_t, double> > Network::neighbors(const size_t &n) const {
	std::vector<std::pair<size_t, double>> neigh;
	for (auto& link : links) {
		if (link.first.first == n) {
			neigh.push_back(std::make_pair(link.first.second, link.second));
		}
	}
	return neigh;
}

std::pair<size_t, double> Network::degree(const size_t &n) const {
	std::pair<size_t, double> sum;
	std::vector<std::pair<size_t, double>> neigh = neighbors(n);
	for (auto& neighbor : neigh) {
		sum.first += 1;
		sum.second += neighbor.second;
	}
	return sum;
}


std::set<size_t> Network::step(const std::vector<double> &thalamic_input) {
	double intensity_excitators, intensity_inhibitors;
	std::set<size_t> firing_neurons;
	double w;
	for (size_t i(0); i < size(); ++ i) {
		std::vector<std::pair<size_t, double>> neighb = neighbors(i);
		for (size_t j(0); j < neighb.size(); ++j) {
			if (neurons[neighb[j].first].firing()) {
				if (neurons[neighb[j].first].is_inhibitory()) {
					intensity_inhibitors += neighb[j].second;
				} else {
					intensity_excitators += neighb[j].second;
				}
			}
		}
		
		
		
		/*std::vector<size_t> neighb;
		for (auto paire : p_neighb) {
			neighb.push_back(paire.first);
		}
		for (auto neigh : neighb) {
			if(neurons[neigh].firing()) {
				if (neurons[neigh].is_inhibitory()) {
					intensity_inhibitors += neighbors(neigh)[i].second;
				} else {
					intensity_excitators += neighbors(neigh)[i].second;
				}
			}
		}*/
		if (neurons[i].is_inhibitory()) {
			w = 0.4;
		} else {
			w = 1;
		}	
		double courant = ((w * thalamic_input[i]) + (0.5 * intensity_excitators) + intensity_inhibitors);
		neurons[i].input(courant);
		if (neurons[i].firing()) {
			firing_neurons.insert(i);
			neurons[i].reset();
		} 
		neurons[i].step();
		intensity_excitators = 0;
		intensity_inhibitors = 0;
	}
	return firing_neurons;
}

