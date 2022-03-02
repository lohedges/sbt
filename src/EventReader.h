/*
  Copyright (c) 2021-2022 Lester Hedges <lester.hedges@gmail.com>

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

#ifndef _EVENTREADER_H
#define _EVENTREADER_H

#include <filesystem>
#include <string>
#include <vector>

struct ModulePair
{
    unsigned hit_start;
    unsigned num_hits;
    float z0;
    float z1;
};

struct Hit
{
    float x;
    float y;
    float z;
    int16_t phi;
    unsigned module_idx;
};

class Event
{
public:
    //! Get the module pairs for this event.
    /*! \return module_pairs
            A vector containing the module pair information for this event.
     */
    const std::vector<ModulePair>& getModulePairs() const;

    //! Get the flattened hits for this event.
    /*! \return hits
            A vector containing the flattened hits for this event, i.e. the
            hits for all module pairs in order.
     */
    const std::vector<Hit>& getHits() const;

    //! Get hits for a specified module pair.
    /*! \param index
            The index of the module pair.

        \return hits
            A vector containing the (sorted by phi) hits for this specified module pair.
     */
    std::vector<Hit> getHits(unsigned index) const;

private:
    friend class EventReader;

    /// The vector of module pairs.
    std::vector<ModulePair> module_pairs;

    /// The vector of hits.
    std::vector<Hit> hits;

    //! Add a module pair.
    /*! \param module_pair
            The module pair record.
     */
    void addModulePair(ModulePair& module_pair);

    //! Add hits for a module pair.
    /*! \param hits
            The hits for a module pair.
     */
    void addHits(std::vector<Hit>& hits);
};

class EventReader
{
public:
    //! Constructor.
    /*! \param path
            The path to the directory containing event data.

        \param verbose
            Whether to print output.
     */
    EventReader(std::filesystem::path path, bool verbose=false);

    //! Read events from file.
    /*! \param verbose
            Whether to print output.
     */
    void readEvents(bool verbose=false);

    //! Get the events.
    /*! \return events
            A vector containing the events.
     */
    const std::vector<Event>& getEvents() const;

private:
    /// A vector of event files.
    std::vector<std::string> file_list;

    /// A vector of events.
    std::vector<Event> events;
};

#endif  /* _EVENTREADER_H */
