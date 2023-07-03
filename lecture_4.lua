#!/usr/bin/env luajit
-- based on https://www.youtube.com/watch?v=y7hJAhKp2d8&feature=youtu.be
require 'ext'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'
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

-- dashboard

-- source parameters
local fmax = 5 * gigahertz	-- maximum frequency

-- grid parameters
local nmax = 1 	-- maximum refractive index.  air = 1.
local NLAM = 20	-- distance between ponits.  "resolve shortest wavelength of highest frequency with 10 points."
				-- "10 is about the threshold where our answers resemble reality."  20 or 40 is better.
local NBUFZ = matrix{100, 100}	-- 100 points before and after the device

-- nominal resolution
local lam0 = c0 / fmax	-- shortest wavelength / free-space wavelength
local dz = lam0/nmax/NLAM	-- grid resolution

-- compute grid size
local Nz = NBUFZ:sum() + 3

-- compute grid axis
local za = matrix(range(0,Nz-1)):emul(dz)

-- set material properties to free space
local ER = matrix{Nz}:ones()
local UR = matrix{Nz}:ones()

-- compute the source

-- compute time step (dt)
local nbc = math.sqrt(UR[1] * ER[1]) -- refractive index at boundaries
		-- sqrt(mu0 * eps0) at the first point
local dt = nbc * dz / (2 * c0)	 -- ensures the wave travels 1 grid cell in 2 time clicks

-- compute source parameters
local tau = .5 / fmax	-- from "theory slides" ... which are where?
local t0 = 5*tau		-- delay

-- compute number of time steps
local tprop = nmax * (Nz * dz) / c0	-- how long it takes the wave to bounce back and forth.
									-- Nz*dz = length of grid
local t = 2 * t0 + 3 * tprop 		-- 2*t0 is the length of the source
									-- 3 * tprop means let the wave bounce 3 times
local STEPS = math.ceil(t/dt)

-- compute the source
local t = matrix(range(0,STEPS-1)) * dt	-- time associated with each timestep
local s = dz / (2 * c0) + dt/2			-- E & M fields are 1/2 timestep and 1/2 grid apart
local nz_src = math.floor(Nz/4 + .5)	-- where to inject the source
local Esrc = (-((t-t0)/tau):epow(2)):map(math.exp)
local A = -math.sqrt(ER[nz_src]/UR[nz_src])
local Hsrc = A*(-((t-t0+s)/tau):epow(2)):map(math.exp)

--[[
gnuplot{
	output = 'tmp.png',
	style = 'data lines',
	data = {t, Esrc},
	{using = '1:2', title='Esrc'},
}
os.exit()
--]]

-- initialize fdtd parameters

-- compute update coefficinets
-- electric field update coefficient  mEy[i] = c0 * dt / ER[i]
local mEy = (matrix{Nz}:ones() * c0 * dt):ediv(ER)
local mHx = (matrix{Nz}:ones() * c0 * dt):ediv(UR)

-- initialize field
local Ey = matrix{Nz}:zeros()
local Hx = matrix{Nz}:zeros()

-- initialize boundary terms
-- these are boundary terms from previous timesteps
local H1 = 0
local H2 = 0
local H3 = 0
local E1 = 0
local E2 = 0
local E3 = 0

local allEys = table()
local allHxs = table()

local loop = (function()
	print'simulating...'
	for T = 1,STEPS do
		-- update H from E
		for nz=1,Nz-1 do
			-- (Hx(i,t+dt) - Hx(i,t))/dt = c0 / ER[i] * (Ey[i+1] - Ey[i]) / dz
			Hx[nz] = Hx[nz] + mHx[nz] * (Ey[nz+1] - Ey[nz]) / dz
		end
		-- right boundary
		Hx[Nz] = Hx[Nz] + mHx[Nz] * (E3 - Ey[Nz]) / dz

		-- H source
		-- correction term pulled out of procedure when doing the total field scattering
		Hx[nz_src-1] = Hx[nz_src-1] - mHx[nz_src-1]*Esrc[T]/dz

		-- record H-field at boundary (of right side)
		H3 = H2
		H2 = H1
		H1 = Hx[1]

		-- update E from H
		Ey[1] = Ey[1] + mEy[1] * (Hx[1] - H3) / dz
		for nz=2,Nz do
			-- (Ey(i,t+dt) - Ey(i,t))/dt = c0 / ER[i] * (Hx[i] - Hx[i-1]) / dz
			Ey[nz] = Ey[nz] + mEy[nz] * (Hx[nz] - Hx[nz-1]) / dz
		end

		-- E source
		Ey[nz_src] = Ey[nz_src] - mEy[nz_src]*Hsrc[T]/dz

		-- record E-field at boundary (of left side)
		E3 = E2
		E2 = E1
		E1 = Ey[Nz]

		-- TODO plot ER, Ey, Hx, dz
		allEys:insert(matrix(Ey))
		allHxs:insert(matrix(Hx))
		coroutine.yield()
	end
	coroutine.yield(false)
end):wrap()

-- [[ using GLApp
local ImGuiApp = require 'imguiapp'
local gl = require 'gl'
local ig = require 'imgui'
local App = class(ImGuiApp)
local max
App.title = 'FTDT test'
function App:update()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	
	max = 5 * math.max((table.sup(Ey)), (table.sup(Hx)))
	gl.glOrtho(-5/Nz, za[Nz]+5/Nz, -max, max, -1, 1)
	
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	local fps = 250
	local thisFrame = math.floor(os.clock() * fps)
	if thisFrame ~= self.lastFrame then
		self.lastFrame = thisFrame
		if loop() == false then loop = function() end end
	end

	for _,info in ipairs{
		{buf=ER, color={1,0,0}},
		{buf=Ey, color={1,1,0}},
		{buf=Hx, color={0,1,1}},
	} do
		gl.glBegin(gl.GL_LINE_STRIP)
		gl.glColor3f(table.unpack(info.color))
		for nz=1,Nz do
			gl.glVertex2f(za[nz], info.buf[nz])
		end
		gl.glEnd()
	end
	App.super.update(self)
end
function App:updateGUI()
	ig.igText('scale '..max)
end
return App():run()
--]]
--[[ gnuplot it all
--assert(dz) print('dz:size()', dz:size())
repeat until not loop()
print'plotting...'
gnuplot{
	persist = true,
	--output = 'out.png',
	style = 'data lines',
	griddata = {x=t, y=za, allEys, allHxs},--, dz},
	--{splot=true, using = '1:2', title='ER'},
	{splot=true, using = '1:2:3', title='Ey'},
	{splot=true, using = '1:2:4', title='Hx'},
	--{using = '1:5', title='dz'},
}
print'done!'
--]]
--[[ interactive 3D plot
repeat until not loop()
local data = table()
for i,ti in ipairs(t) do
	for j,zj in ipairs(za) do
		data:insert{ti, zj, allEys[i][j], allHxs[i][j]}
	end
end
local t, z, Ey, Hx = matrix(data):transpose():unpack()
require 'plot3d'{
	Ey = {t, z, Ey, enabled=true}, 
	Hx = {t, z, Hx, enabled=true},
}
--]]
