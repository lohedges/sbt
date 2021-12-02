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

#include <iostream>

#include <poplar/DeviceManager.hpp>
#include <poplar/IPUModel.hpp>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

/// Setup and connect to a Poplar IPU device.
poplar::Device setIpuDeviceWithIpuModelFallBack(bool &is_model);

int main()
{
    // Read events from file.
    EventReader event_reader(std::filesystem::path("data"), true);
    event_reader.readEvents(true);

    // Setup the IPU device handle.
    bool is_model;
    auto device = setIpuDeviceWithIpuModelFallBack(is_model);

    // Setup the graph program.
    std::cout << "Creating graph program...\n";
    SearchByTripletIPU search_by_triplet(std::move(device),
                                         event_reader.getEvents(),
                                         is_model);

    double events_per_sec;
    std::cout << "Running benchmarks...\n";
    // Run the program.
    std::tuple<unsigned*, Track*, unsigned*, Tracklet*> result =
        search_by_triplet.execute(events_per_sec);

    // Output timing statistics.
    std::cout << events_per_sec << " event/s" << std::endl;

    return 0;
}

poplar::Device setIpuDeviceWithIpuModelFallBack(bool &is_model)
{
    auto dm = poplar::DeviceManager::createDeviceManager();
    auto hwDevices = dm.getDevices(poplar::TargetType::IPU, 1);
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
