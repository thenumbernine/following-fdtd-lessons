#!/usr/bin/env luajit
-- based on https://www.youtube.com/watch?v=y7hJAhKp2d8&feature=youtu.be
local matrix = require 'matrix'

local meters = 1
local centimeters = 1e-2 * meters
local millimeters = 1e-3 * meters
local inches = 2.54 * centimeters
local feet = 12 * inches
local seconds = 1
local hertz = 1/seconds
local kilohertz = 1e+3 * hertz
local megahertz = 1e+6 * hertz
local gigahertz = 1e+9 * hertz
-- hmm, should henries and farads both be 1?
local farads = 1
local henries = 1
local c0 = 299792458 * meters / seconds
local e0 = 8.8541878176e-12 * farads/meters
local u0 = 1.2566370614e-6 * henries/meters

local dz = .006 * meters
local Nz = 200
local dt = 1e-11 * seconds
local STEPS = 1000

-- set material properties to free space
local ER = matrix{Nz}:ones()
local UR = matrix{Nz}:ones()

-- electric field update coefficient  mEy[i] = c0 * dt / ER[i]
local mEy = (matrix{Nz}:ones() * c0 * dt):ediv(ER)
local mHx = (matrix{Nz}:ones() * c0 * dt):ediv(UR)

-- field
local Ey = matrix{Nz}:zeros()
local Hx = matrix{Nz}:zeros()

for T = 1,STEPS do
	-- update H from E
	for nz=1,Nz-1 do
		-- (Hx(i,t+dt) - Hx(i,t))/dt = c0 / ER[i] * (Ey[i+1] - Ey[i]) / dz
		Hx[nz] = Hx[nz] + mHx[nz] * (Ey[nz+1] - Ey[nz]) / dz
	end
	-- right boundary
	local Ey_rhs = 0
	Hx[Nz] = Hx[Nz] + mHx[Nz] * (Ey_rhs - Ey[Nz]) / dz

	-- update E from H
	local Hx_lhs = 0
	Ey[1] = Ey[1] + mEy[1] * (Hx[1] - Hx_lhs) / dz
	for nz=2,Nz do
		-- (Ey(i,t+dt) - Ey(i,t))/dt = c0 / ER[i] * (Hx[i] - Hx[i-1]) / dz
		Ey[nz] = Ey[nz] + mEy[nz] * (Hx[nz] - Hx[nz-1]) / dz
	end

	-- TODO plot ER, Ey, Hz, dz
end
