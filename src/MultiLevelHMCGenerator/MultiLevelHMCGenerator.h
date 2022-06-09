/**
 * @file       MultiLevelHMCGenerator.h
 * @brief      Declarations of Multi level HMC
 * @author     nico
 * @version    0.0.2
 * @date       26.03.22
 */


#ifndef BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
#define BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H

#include <BaseModel.h>
#include <HMCGenerator.h>
#include <Analyzer.h>
#include <random>
#include <iostream>
#include <chrono>

#include <utility>
#include <memory>

/**
 * @brief Template of the Multilevel HMC algorithm
 * @tparam configuration_type Datatype, which is used for the configurations of the model
 */
template<class configuration_type>
class MultiLevelHMCGenerator {
public:
    /**
     * @brief Constructor of the Multi Level HMC generator
     * @param model_ Model for which to generate ensembles
     * @param nu_pre_ Amount of pre coarsening steps to take at each level
     * @param nu_post_ Amount of post coarsening steps to take at each level
     * @param gamma_ Amount of repetitions at each level (determines if a 'V' or 'W' or ... cycle is performed)
     * @param InterpolationType_ Interpolation type used to generate the coarser levels
     * @param amount_of_steps_ Amount of steps to be used in the integration process for each level
     * @param step_sizes_ Step size in the integration process for each level
     * @param generator_ Random number generator to be used for the HMC process
     */
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, std::vector<int> erg_jump_dists_, size_t gamma_,
                           InterpolationType InterpolationType_,
                           const std::vector<size_t> &amount_of_steps_, const std::vector<double> &step_sizes_,
                           std::default_random_engine &generator_);

    /**
     * @brief Loads the MultiLevelHMCGenerator from \ p file
     * @param model_ Model for which to generate ensembles
     * @param file File which has the levels as groups, which hold the level attributes
     * @param generator_ Random number generator to be used for the HMC process
     */
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                           HighFive::File &file,
                           std::default_random_engine &generator_);

    /**
     * @brief Generate \p amount_of_samples amount of ensembles, starting from \p phiStart and doing
     *        \p amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rates
     */
    std::vector<double> generate_ensembles(const configuration_type &phiStart,
                                           size_t amount_of_samples, size_t amount_of_thermalization_steps = 10);

    /**
     * @brief Calculates and dumps the observable \a observable_function_pointer for each level to \p file
     *        into subgroups for each level intern into the \p name dataset
     * @param observable_function_pointer Pointer to function, which returns the value of the observable
     *                                    for the given configuration
     * @param name Name under which the dataset will be stored
     * @param file File to dump to
     */
    void dump_observable(double (
    BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
                         const std::string &name, HighFive::File &file);

    void analyze_dataset(const std::string &name, HighFive::File &file,
                         int block_size, int size_to_use, size_t start_index, int amount_of_sample_sets, size_t max_t);

    /**
     * @brief Dumps the MultiLevelHMCGenerator, including all sub models and sub HMCGenerators
     *        into different subgroups with corresponding level parameters
     * @param file File to dump to
     */
    void dumpToH5(HighFive::File &file);

    /**
     * @brief Propagates the change of attributes from the fine level model to all coarser levels
     */
    [[maybe_unused]] void propagate_update();


private:
    /**
     * @brief Recursion to go to coarser levels
     * @param level Current level id (finest=0, coarsest=\c nu_pre.size()-1 )
     * @param phi Starting configuration/field
     * @return Updated configuration/field
     */
    configuration_type LevelRecursion(int level, const configuration_type &phi);

    /**
     * @brief Prefix for the level groups in H5 files
     */
    [[maybe_unused]] static const char *level_name;

    /**
     * @brief Name for the measurements groups in H5 files
     */
    [[maybe_unused]] static const char *measurements_name;

    /**
     * @brief Name for the data dataset of measurements in H5 files
     */
    [[maybe_unused]] static const char *measurements_data_name;

    /**
     * @brief Amount of pre coarsening steps to take at each level
     */
    std::vector<size_t> nu_pre;
    /**
     * @brief String to be used as key for the \a nu_pre in H5 files
     */
    [[maybe_unused]] static const char *nu_pre_name;

    /**
     * @brief Amount of post coarsening steps to take at each level
     */
    std::vector<size_t> nu_post;
    /**
     * @brief String to be used as key for the \a nu_post in H5 files
     */
    [[maybe_unused]] static const char *nu_post_name;

    /**
     * @brief Distance between ergodicity jumps at each level
     */
    std::vector<int> erg_jump_dists;
    /**
     * @brief String to be used as key for the \a erg_jump_dists in H5 files
     */
    [[maybe_unused]] static const char *erg_jump_dists_name;

    /**
     * @brief Amount of repetitions at each level (determines if a 'V' or 'W' or ... cycle is performed)
     */
    size_t gamma;
    /**
     * @brief String to be used as key for the \a gamma in H5 files
     */
    [[maybe_unused]] static const char *gamma_name;

    /**
     * @brief Amount of ticks in milliseconds taken for the production
     */
    size_t tick_time;
    /**
     * @brief String to be used as key for the \a tick_time in H5 files
     */
    [[maybe_unused]] static const char *tick_time_name;

    /**
     * @brief Interpolation type used to generate the coarser levels
     */
    InterpolationType inter_type;
    /**
     * @brief String to be used as key for the \a inter_type in H5 files
     */
    [[maybe_unused]] static const char *inter_type_name;

    /**
     * @brief Random number generator to be used for the HMC process
     */
    std::default_random_engine &generator;

    /**
     * @brief Stores the HMC generator for each level
     */
    std::vector<HMCGenerator<configuration_type>> HMCStack;

    /**
     * @brief Stores the model for each level
     */
    std::vector<std::unique_ptr<BaseModel<configuration_type>>> ModelStack;

    /**
     * @brief Saves the Acceptance rate at each level
     */
    std::vector<double> AcceptanceRates;
    /**
     * @brief String to be used as key for the \a AcceptanceRate in H5 files
     */
    [[maybe_unused]] static const char *AcceptanceRate_name;
};

template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::level_name{"level"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        measurements_name{"measurements"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        measurements_data_name{"data"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::nu_pre_name{"nu_pre"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::nu_post_name{"nu_post"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        erg_jump_dists_name{"erg_jump_dists"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::gamma_name{"gamma"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::tick_time_name{"tick_time"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        inter_type_name{"inter_type"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        AcceptanceRate_name{"AcceptanceRate"};


template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                                                                   std::vector<size_t> nu_pre_,
                                                                   std::vector<size_t> nu_post_,
                                                                   std::vector<int> erg_jump_dists_,
                                                                   size_t gamma_,
                                                                   InterpolationType InterpolationType_,
                                                                   const std::vector<size_t> &amount_of_steps_,
                                                                   const std::vector<double> &step_sizes_,
                                                                   std::default_random_engine &generator_)
        : nu_pre{std::move(nu_pre_)}, nu_post{std::move(nu_post_)}, erg_jump_dists{std::move(erg_jump_dists_)},
          gamma{gamma_}, tick_time{}, inter_type{InterpolationType_}, generator{generator_}, AcceptanceRates{} {
    //TODO: add auto sizing
    bool size_ok = (gamma > 0)
                   && (nu_pre.size() == nu_post.size())
                   && (nu_pre.size() == erg_jump_dists.size())
                   && (nu_pre.size() == amount_of_steps_.size())
                   && (nu_pre.size() == step_sizes_.size());
    assert(size_ok);
    if (!size_ok) {
        std::cerr << "Problem with sizes!\n";
        exit(-42);
    }


    AcceptanceRates.resize(nu_pre.size());
    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model_.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(*ModelStack[0], amount_of_steps_[0], step_sizes_[0], generator));

    assert(nu_post[0] > 0);
    for (int i = 1; i < nu_pre.size(); ++i) {
        assert(nu_post[i] > 0);
        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_coarser_model(inter_type)));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], amount_of_steps_[i], step_sizes_[i], generator));
    }

}

template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                                                                   HighFive::File &file,
                                                                   std::default_random_engine &generator_)
        : gamma{}, tick_time{}, generator{generator_}, AcceptanceRates{} {

    std::string current_level_name{level_name};
    current_level_name.append(std::to_string(0));
    assert(file.exist(current_level_name));

    HighFive::Group current_level = file.getGroup(current_level_name);

    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model_.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(*ModelStack[0], current_level, generator));

    size_t buffer;
    current_level.getAttribute(nu_pre_name).read(buffer);
    nu_pre.push_back(buffer);
    current_level.getAttribute(nu_post_name).read(buffer);
    nu_post.push_back(buffer);
    int iBuffer;
    if (current_level.hasAttribute(erg_jump_dists_name)) {
        current_level.getAttribute(erg_jump_dists_name).read(iBuffer);
    } else {
        iBuffer = -1;
    }
    erg_jump_dists.push_back(iBuffer);
    current_level.getAttribute(gamma_name).read(gamma);
    if (current_level.hasAttribute(tick_time_name)) {
        current_level.getAttribute(tick_time_name).read(tick_time);
    }
    current_level.getAttribute(inter_type_name).read(iBuffer);
    inter_type = static_cast<InterpolationType>(iBuffer);
    double dBuffer;
    current_level.getAttribute(AcceptanceRate_name).read(dBuffer);
    AcceptanceRates.push_back(dBuffer);
    assert(nu_post[0] > 0);

    for (int i = 1; i < 1000; ++i) {
        current_level_name = level_name;
        current_level_name.append(std::to_string(i));
        if (!file.exist(current_level_name)) {
            break;
        }

        current_level = file.getGroup(current_level_name);


        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_model_at(current_level)));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], current_level, generator));

        current_level.getAttribute(nu_pre_name).read(buffer);
        nu_pre.push_back(buffer);
        current_level.getAttribute(nu_post_name).read(buffer);
        nu_post.push_back(buffer);
        if (current_level.hasAttribute(erg_jump_dists_name)) {
            current_level.getAttribute(erg_jump_dists_name).read(iBuffer);
        } else {
            iBuffer = -1;
        }
        erg_jump_dists.push_back(iBuffer);
        current_level.getAttribute(gamma_name).read(gamma);
        if (current_level.hasAttribute(tick_time_name)) {
            current_level.getAttribute(tick_time_name).read(tick_time);
        }
        current_level.getAttribute(inter_type_name).read(iBuffer);
        inter_type = static_cast<InterpolationType>(iBuffer);
        current_level.getAttribute(AcceptanceRate_name).read(dBuffer);
        AcceptanceRates.push_back(dBuffer);
        assert(nu_post[i] > 0);
    }
}

template<class configuration_type>
std::vector<double> MultiLevelHMCGenerator<configuration_type>::generate_ensembles(const configuration_type &phiStart,
                                                                                   size_t amount_of_samples,
                                                                                   size_t amount_of_thermalization_steps) {
    configuration_type phi(phiStart);
    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        if (i % (amount_of_thermalization_steps / 10) == 0) {
            std::cout << "|";
            std::cout.flush();
        }
        if (i % (amount_of_thermalization_steps / 100) == 0) {
            std::cout << "=";
            std::cout.flush();
        }
        phi = LevelRecursion(0, phi);
    }
    std::cout << std::endl;
    HMCStack[0].clear_ensembles();
    std::cout << "Thermalization done." << std::endl;
    for (int i = 0; i < AcceptanceRates.size(); ++i) {
        std::cout << "Acceptance rate:" << AcceptanceRates[i] << "\t\t" <<
                  AcceptanceRates[i] / (amount_of_thermalization_steps * (nu_pre[i] + nu_post[i]) * int_pow(gamma, i))
                  << std::endl;
        AcceptanceRates[i] = 0.;
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < amount_of_samples; ++i) {
        if (i % (amount_of_samples / 10) == 0) {
            std::cout << "|";
            std::cout.flush();
        }
        if (i % (amount_of_samples / 100) == 0) {
            std::cout << "=";
            std::cout.flush();
        }
        phi = LevelRecursion(0, phi);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    tick_time = duration.count();

    std::cout << std::endl;
    std::cout << "Production done.\nExecution time: " << tick_time / 1e3 << " s" << std::endl;

    for (int i = 0; i < AcceptanceRates.size(); ++i) {
        AcceptanceRates[i] = AcceptanceRates[i] / (amount_of_samples * (nu_pre[i] + nu_post[i]) * int_pow(gamma, i));
    }
    return AcceptanceRates;
}

template<class configuration_type>
configuration_type
MultiLevelHMCGenerator<configuration_type>::LevelRecursion(int level, const configuration_type &phi) {
    configuration_type currentField{phi};
    AcceptanceRates[level] += HMCStack[level].generate_ensembles(currentField, nu_pre[level], 0,
                                                                 level == 0, erg_jump_dists[level]);
    currentField = HMCStack[level].get_last_configuration(currentField);
    //std::cout << "Level:\t" << level << std::endl;

    if (level < (nu_pre.size() - 1)) {
        ModelStack[level + 1]->update_fields(currentField);
        configuration_type CoarseCorrections = ModelStack[level + 1]->get_empty_field();
        for (int i = 0; i < gamma; ++i) {
            CoarseCorrections = LevelRecursion(level + 1, CoarseCorrections);
        }
        ModelStack[level + 1]->interpolate(CoarseCorrections, currentField);
    }
    AcceptanceRates[level] += HMCStack[level].generate_ensembles(currentField, nu_post[level], 0,
                                                                 level == 0, erg_jump_dists[level]);
    return HMCStack[level].get_last_configuration(currentField);
}

template<class configuration_type>
[[maybe_unused]] void MultiLevelHMCGenerator<configuration_type>::propagate_update() {
    //Start from 0 because this way also level 0 updates its parameters from the root model
    for (int i = 0; i < ModelStack.size(); ++i) {
        ModelStack[i]->pull_attributes_from_root();
    }
}


template<class configuration_type>
void MultiLevelHMCGenerator<configuration_type>::dumpToH5(HighFive::File &file) {
    for (int i = 0; i < HMCStack.size(); ++i) {
        HighFive::Group current_level = file.getGroup(file.getPath());

        std::string current_level_name{level_name};
        current_level_name.append(std::to_string(i));

        if (file.exist(current_level_name)) {
            current_level = file.getGroup(current_level_name);
        } else {
            current_level = file.createGroup(current_level_name);
        }
        HMCStack[i].dumpToH5(current_level);
        write_static_size(nu_pre[i], current_level, nu_pre_name);
        write_static_size(nu_post[i], current_level, nu_post_name);
        write_static_size(erg_jump_dists[i], current_level, erg_jump_dists_name);
        write_static_size(gamma, current_level, gamma_name);
        write_static_size(tick_time, current_level, tick_time_name);
        write_static_size(static_cast<int>(inter_type), current_level, inter_type_name);
        write_static_size(AcceptanceRates[i], current_level, AcceptanceRate_name);
    }
}

template<class configuration_type>
void MultiLevelHMCGenerator<configuration_type>::dump_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
        const std::string &name, HighFive::File &file) {
    for (int i = 0; i < HMCStack.size(); ++i) {
        HighFive::Group current_level = file.getGroup(file.getPath());
        HighFive::Group measurements = current_level;
        HighFive::Group observable_group = current_level;

        std::string current_level_name{level_name};
        current_level_name.append(std::to_string(i));

        if (file.exist(current_level_name)) {
            current_level = file.getGroup(current_level_name);
        } else {
            current_level = file.createGroup(current_level_name);
        }
        if (current_level.exist(measurements_name)) {
            measurements = current_level.getGroup(measurements_name);
        } else {
            measurements = current_level.createGroup(measurements_name);
        }
        if (measurements.exist(name)) {
            observable_group = measurements.getGroup(name);
        } else {
            observable_group = measurements.createGroup(name);
        }
        HMCStack[i].dumpToH5(current_level, true);
        write_static_size(nu_pre[i], current_level, nu_pre_name);
        write_static_size(nu_post[i], current_level, nu_post_name);
        write_static_size(erg_jump_dists[i], current_level, erg_jump_dists_name);
        write_static_size(gamma, current_level, gamma_name);
        write_static_size(tick_time, current_level, tick_time_name);
        write_static_size(static_cast<int>(inter_type), current_level, inter_type_name);
        write_static_size(AcceptanceRates[i], current_level, AcceptanceRate_name);
        HighFive::DataSet dataset = HMCStack[i].dump_observable(observable_function_pointer, measurements_data_name,
                                                                observable_group);
    }
}

template<class configuration_type>
void MultiLevelHMCGenerator<configuration_type>::analyze_dataset(const std::string &name, HighFive::File &file,
                                                                 int block_size, int size_to_use, size_t start_index,
                                                                 int amount_of_sample_sets,
                                                                 size_t max_t) {
    for (int i = 0; i < HMCStack.size(); ++i) {
        std::string current_level_name{level_name};
        current_level_name.append(std::to_string(i));
        if (!file.exist(current_level_name)) {
            std::cerr << "No group named " << current_level_name << '\n';
            return;
        }
        HighFive::Group current_level = file.getGroup(current_level_name);

        if (!current_level.exist(measurements_name)) {
            std::cerr << "No group named " << measurements_name << '\n';
            return;
        }
        HighFive::Group measurements = current_level.getGroup(measurements_name);

        if (!measurements.exist(name)) {
            std::cerr << "No group named " << name << '\n';
            return;
        }
        HighFive::Group observable_group = measurements.getGroup(name);

        if (!observable_group.exist(measurements_data_name)) {
            std::cerr << "No dataset named " << measurements_data_name << '\n';
            return;
        }
        if (i == 0) {
            Analyzer a(observable_group, measurements_data_name, start_index, generator);
            a.auto_correlation(max_t);
            a.block_data(block_size, size_to_use);
            a.bootstrap_data(amount_of_sample_sets);
        }
    }
}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
