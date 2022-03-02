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

#include <iostream>
#include <vector>

#include <poplar/DeviceManager.hpp>
#include <poplar/IPUModel.hpp>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

/// Setup and connect to a Poplar IPU device.
poplar::Device setIpuDeviceWithIpuModelFallBack(bool &is_model, const int num_ipus);

int main(int argc, char *argv[])
{
    // Set the default number of IPUs to run on.
    int num_ipus = 1;

    // Get the number of IPUs from the command-line.
    if (argc > 1)
    {
        std::string arg = argv[1];

        try
        {
            std::size_t pos;
            num_ipus = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number of IPUs: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number of IPUs: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number of IPUs out of range: " << arg << '\n';
            exit(-1);
        }
    }

    if (num_ipus > 4)
    {
        std::cerr << "Number of IPUs must be in range 1-4.\n";
        exit(-1);
    }

    // Read events from file.
    EventReader event_reader(std::filesystem::path("data"), true);
    event_reader.readEvents(true);

    // Setup the IPU device handle.
    bool is_model;
    auto device = setIpuDeviceWithIpuModelFallBack(is_model, num_ipus);

    if (is_model and num_ipus > 1)
    {
        std::cout << "Running on an IPUModel. Ignoring 'num_ipus' option!\n";
    }

    // Setup the graph program.
    std::cout << "Creating graph program...\n";
    SearchByTripletIPU search_by_triplet(std::move(device),
                                         event_reader.getEvents(),
                                         is_model);

    if (num_ipus == 1)
    {
        std::cout << "Running benchmarks on 1 IPU...\n";
    }
    else
    {
        std::cout << "Running benchmarks on " << num_ipus << "IPUs...\n";
    }

    // Run the program.
    std::tuple<std::vector<unsigned>, std::vector<Track>,
               std::vector<unsigned>, std::vector<Tracklet>>
        result = search_by_triplet.execute();

    return 0;
}

poplar::Device setIpuDeviceWithIpuModelFallBack(bool &is_model, const int num_ipus)
{
    auto dm = poplar::DeviceManager::createDeviceManager();
    auto hwDevices = dm.getDevices(poplar::TargetType::IPU, num_ipus);
    if (hwDevices.size() > 0)
    {
        for (auto &d : hwDevices)
        {
            if (d.attach())
            {
                is_model = false;
                return std::move(d);
            }
        }
    }

    poplar::IPUModel ipuModel;

    is_model = true;
    return ipuModel.createDevice();
}
