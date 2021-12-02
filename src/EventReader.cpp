/*
  Copyright (c) 2021 Lester Hedges <lester.hedges@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <iostream>
#include <sstream>

#include "Definitions.h"
#include "EventReader.h"

EventReader::EventReader(std::filesystem::path path, bool verbose)
{
    // Make sure the event data directory is valid.
    if(not std::filesystem::is_directory(path))
    {
        throw std::runtime_error("Invalid directory: " + path.string());
    }

    // Iterate over all files in the directory and store their names.
    if (verbose)
        std::cout << "Scanning event directory..." << std::endl;
    for (const auto &file : std::filesystem::directory_iterator(path))
    {
        const auto full_name = file.path().string();

        if (file.is_regular_file())
        {
            this->file_list.push_back(full_name);
            if (verbose)
                std::cout << "  " << full_name << std::endl;
        }
        else
        {
            throw std::runtime_error("Invalid event file: " + full_name);
        }
    }
}

void EventReader::readEvents(bool verbose)
{
    // Loop over all event files and parse into an Event object.
    if (verbose)
        std::cout << "Reading event files..." << std::endl;
    for (const auto &file : this->file_list)
    {
        if (verbose)
            std::cout << "  " << file << std::endl;

        std::ifstream event_file(file);

        if (event_file.is_open())
        {
            bool is_hit = false;
            std::string line;

            Event event;
            ModulePair module_pair;
            Hit hit;
            std::vector<Hit> hits;

            // Index of the starting hit in the flattened hits array for the
            // current module pair.
            unsigned hit_start = 0;

            // The number of module pairs recorded for this event.
            unsigned num_mp = 0;

            while (std::getline(event_file, line))
            {
                // End of module pair record.
                if (line.empty())
                {
                    // Store the information.
                    event.addModulePair(module_pair);
                    event.addHits(hits);

                    num_mp++;
                    hit_start += module_pair.num_hits;
                    is_hit = false;
                    hits.clear();

                    continue;
                }

                std::stringstream ss(line);

                // Store the module pair record.
                if (not is_hit)
                {
                    module_pair.hit_start = hit_start;

                    ss >> module_pair.num_hits
                       >> module_pair.z0
                       >> module_pair.z1;

                    is_hit = true;
                    continue;
                }

                // Now process the hits for this module pair.
                if (is_hit)
                {
                    ss >> hit.x
                       >> hit.y
                       >> hit.z
                       >> hit.phi
                       >> hit.module_idx;

                    hits.push_back(hit);
                }
            }

            event_file.close();

            if (num_mp != constants::num_module_pairs)
            {
                std::string msg = "Event has incorrect number of module pairs. ";
                msg += "Found " + std::to_string(num_mp) + ", expected "
                                + std::to_string(constants::num_module_pairs);
                throw std::runtime_error(msg);
            }

            // Store the event.
            this->events.push_back(event);
        }
        else
        {
            throw std::runtime_error("Unable to read event file: " + file);
        }
    }
}

const std::vector<Event>& EventReader::getEvents() const
{
    return this->events;
}

void Event::addModulePair(ModulePair& module_pair)
{
    this->module_pairs.push_back(module_pair);
}

void Event::addHits(std::vector<Hit>& hits)
{
    this->hits.insert(this->hits.end(), hits.begin(), hits.end());
}

const std::vector<ModulePair>& Event::getModulePairs() const
{
    return this->module_pairs;
}

const std::vector<Hit>& Event::getHits() const
{
    return this->hits;
}

std::vector<Hit> Event::getHits(unsigned index) const
{
    if (index < 0 or index >= this->module_pairs.size())
    {
        throw std::runtime_error("Module pair index out of bounds!");
    }

    // Work out the number of hits up to this module pair.
    unsigned int num_hits = 0;
    for (unsigned i=0; i<index; ++i)
        num_hits += this->module_pairs[i].num_hits;

    // Extract the hits for this module pair from the flattened vector.
    std::vector<Hit>::const_iterator first
        = this->hits.begin() + num_hits;
    std::vector<Hit>::const_iterator last
        = this->hits.begin() + num_hits + this->module_pairs[index].num_hits;
    std::vector<Hit> module_pair_hits(first, last);

    return module_pair_hits;
}
