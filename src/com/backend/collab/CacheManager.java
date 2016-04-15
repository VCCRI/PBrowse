package com.backend.collab;

import net.sf.ehcache.Cache;
import net.sf.ehcache.config.CacheConfiguration;

/**
 * Publicly exposes the cache manager so we can invalidate entries on relevant database changes
 * @author root
 *
 */
public class CacheManager 
{
	public Cache ddCache 			= null;
	public Cache groupCache 		= null;
	public Cache userfilesCache		= null;
	public Cache userCommentsCache	= null;
    
	private static CacheManager cm 	= null;
	
	private CacheManager()
	{
		//initialize the cache manager for file DD requests
		net.sf.ehcache.CacheManager singletonManager = net.sf.ehcache.CacheManager.create();
        
		//init 2 caches - DD cache stores DDs of requested files
		//group cache stores user membership entries
		//both are used exclusively by the file servlet
		singletonManager.addCache("ddCache");
		singletonManager.addCache("groupCache");
		singletonManager.addCache("userfilesCache");
		singletonManager.addCache("userCommentsCache");
        
		ddCache = singletonManager.getCache("ddCache");
		groupCache = singletonManager.getCache("groupCache");
		userfilesCache = singletonManager.getCache("userfilesCache");
		userCommentsCache = singletonManager.getCache("userCommentsCache");
        
		//ensure cache never times out cache elements are only invalidated by certain DBA calls
		//which modify the underlying structures
		CacheConfiguration config = ddCache.getCacheConfiguration();
        config.setTimeToIdleSeconds(Long.MAX_VALUE);
        config.setTimeToLiveSeconds(Long.MAX_VALUE);
        
        config = groupCache.getCacheConfiguration();
        config.setTimeToIdleSeconds(Long.MAX_VALUE);
        config.setTimeToLiveSeconds(Long.MAX_VALUE);
        
        config = userfilesCache.getCacheConfiguration();
        config.setTimeToIdleSeconds(Long.MAX_VALUE);
        config.setTimeToLiveSeconds(Long.MAX_VALUE);
        
        config = userCommentsCache.getCacheConfiguration();
        config.setTimeToIdleSeconds(Long.MAX_VALUE);
        config.setTimeToLiveSeconds(Long.MAX_VALUE);
	}
	
	public static CacheManager cm()
	{
		if (cm == null)
		{
			cm = new CacheManager();
		}
		return cm;
	}
	
}
